import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid, SpinLock

    #  ----------------------------------------------------------------------------
    import Chemostat_Heerden2013
    const ChH = Chemostat_Heerden2013

    const iJR = ChH.iJR904
    const Hd = ChH.HeerdenData # experimental data
    const Bd = ChH.BegData    # cost data

    #  ----------------------------------------------------------------------------
    # run add "https://github.com/josePereiro/Chemostat" in the 
    # julia Pkg REPL for installing the package
    import Chemostat
    import Chemostat.LP.MathProgBase
    const Ch = Chemostat
    const ChU = Ch.Utils
    const ChLP = Ch.LP

    using Serialization

    import UtilsJL
    const UJL = UtilsJL

    using Base.Threads
end

## ----------------------------------------------------------------------------------
DAT = ChU.DictTree();

## ----------------------------------------------------------------------------------
FLX_IDERS = ["GLC", "SUCC", "AC", "FORM"]
DAT[:FLX_IDERS] = FLX_IDERS;

## -------------------------------------------------------------------
# ME methods
const ME_Z_OPEN_G_OPEN        = :ME_Z_OPEN_G_OPEN
const ME_MAX_POL              = :ME_MAX_POL
const ME_Z_EXPECTED_G_MOVING  = :ME_Z_EXPECTED_G_MOVING
const ME_Z_EXPECTED_G_BOUNDED = :ME_Z_EXPECTED_G_BOUNDED
const ME_Z_FIXXED_G_BOUNDED   = :ME_Z_FIXXED_G_BOUNDED

# LP methods
const FBA_Z_FIX_MIN_COST      = :FBA_Z_FIX_MIN_COST
const FBA_MAX_BIOM_MIN_COST   = :FBA_MAX_BIOM_MIN_COST
const FBA_Z_FIX_MIN_VG_COST   = :FBA_Z_FIX_MIN_VG_COST
const FBA_Z_VG_FIX_MIN_COST   = :FBA_Z_VG_FIX_MIN_COST

LP_METHODS = [
    FBA_Z_FIX_MIN_COST, FBA_MAX_BIOM_MIN_COST, 
    FBA_Z_FIX_MIN_VG_COST, FBA_Z_VG_FIX_MIN_COST
]
DAT[:LP_METHODS] = LP_METHODS

ME_METHODS = [
    # ME_Z_OPEN_G_OPEN, 
    ME_MAX_POL,
    # ME_Z_FIXXED_G_BOUNDED, 
    # ME_Z_EXPECTED_G_BOUNDED, 
    # ME_Z_EXPECTED_G_MOVING
]
DAT[:ME_METHODS] = ME_METHODS

ALL_METHODS = [LP_METHODS; ME_METHODS]
DAT[:ALL_METHODS] = ALL_METHODS

EXPS = Hd.EXPS 
DAT[:EXPS] = EXPS

## -------------------------------------------------------------------
ME_INDEX_FILE = iJR.procdir("maxent_ep_index.bson")
ME_INDEX = ChU.load_data(ME_INDEX_FILE; verbose = false);

## -------------------------------------------------------------------
LP_DAT_FILE = iJR.procdir("lp_dat_file.bson")
LP_DAT = ChU.load_data(LP_DAT_FILE; verbose = false);

## ----------------------------------------------------------------------------------
Hd_mets_map = iJR.load_mets_map()
Hd_rxns_map = iJR.load_rxns_map()

## ----------------------------------------------------------------------------------
# COMMON DAT
let
    max_model = iJR.load_model("max_model"; uncompress = false)
    objider = iJR.BIOMASS_IDER

    for exp in EXPS
        # exp dat
        Hd_biom = Hd.val("D", exp)
        DAT[:exp, :flx, "D", exp] = Hd_biom
        DAT[:exp, :err, "D", exp] = 0.0

        # bounds
        fva_model = iJR.load_model("fva_models", exp; uncompress = false)
        max_lb, max_ub = ChU.bounds(max_model, objider)
        fva_lb, fva_ub = ChU.bounds(fva_model, objider)
        lb = max(max_lb, fva_lb)
        ub = min(max_ub, fva_ub)
        DAT[:bounds, "D", exp] = (lb, ub)

        for Hd_ider in FLX_IDERS
            model_exch = Hd_rxns_map[Hd_ider]

            # exp dat
            Hd_flx = Hd.uval(Hd_ider, exp)
            DAT[:exp, :flx, Hd_ider, exp] = Hd_flx
            DAT[:exp, :err, Hd_ider, exp] = 0.0

            # bounds
            max_lb, max_ub = ChU.bounds(max_model, model_exch)
            fva_lb, fva_ub = ChU.bounds(fva_model, model_exch)
            lb = max(max_lb, fva_lb)
            ub = min(max_ub, fva_ub)
            DAT[:bounds, Hd_ider, exp] = (lb, ub)

        end
    end
end

## ----------------------------------------------------------------------------------
# MAXENT DAT
let 
    WLOCK = ReentrantLock()
    objider = iJR.BIOMASS_IDER

    # Feed jobs
    nths = nthreads()
    Ch = Channel(nths) do ch
        for exp in DAT[:EXPS], method in DAT[:ME_METHODS]
            put!(ch, (exp, method))
        end
    end

    @threads for thid in 1:nths
        thid = threadid()
        thid == 1 && nths > 1 && continue
        for (exp, method) in Ch
            
            # ME data
            datfile = ME_INDEX[method, :DFILE, exp]
            datfile == :unfeasible && continue
            dat = deserialize(datfile)
            
            model = dat[:model]
            epouts = dat[:epouts]
            exp_beta = maximum(keys(epouts)) # dat[:exp_beta]
            epout = epouts[exp_beta]
            
            lock(WLOCK) do
                @info("Doing", 
                    exp, method, 
                    length(dat[:epouts]), 
                    epout.iter, thid
                ); println()
            end

            # Biomass
            ep_biom = ChU.av(model, epout, objider)
            ep_std = sqrt(ChU.va(model, epout, objider))
            
            # store
            lock(WLOCK) do
                DAT[method, :flx, "D", exp] = ep_biom
                DAT[method, :err, "D", exp] = ep_std
            end

            for Hd_met in FLX_IDERS
                model_exch = Hd_rxns_map[Hd_met]

                # flxs
                ep_av = ChU.av(model, epout, model_exch)
                ep_std = sqrt(ChU.va(model, epout, model_exch))

                # proj 2d
                proj = ChLP.projection2D(model, objider, model_exch; l = 50)
                        
                lock(WLOCK) do
                    DAT[method, :proj, Hd_met, exp] = proj
                    DAT[method, :flx, Hd_met, exp] = ep_av
                    DAT[method, :err, Hd_met, exp] = ep_std
                end
            end

        end # for (exp, method)
    end # for thid
end

## ----------------------------------------------------------------------------------
# LP DAT
let
    objider = iJR.BIOMASS_IDER

    for method in LP_METHODS
            
        for exp in EXPS

            model = LP_DAT[method, :model, exp]
            fbaout = LP_DAT[method, :fbaout, exp]

            # Biomass
            fba_flx = ChU.av(model, fbaout, objider)
            Hd_flx = Hd.val("D", exp)
            DAT[method, :flx, "D", exp] = fba_flx

            for Hd_ider in FLX_IDERS
                model_ider = Hd_rxns_map[Hd_ider]

                Hd_flx = Hd.uval(Hd_ider, exp)
                fba_flx = ChU.av(model, fbaout, model_ider)
                DAT[method, :flx, Hd_ider, exp] = fba_flx
            end
        end

    end # for method
end

## ----------------------------------------------------------------------------------
DAT_FILE = iJR.procdir("dat.bson")
UJL.save_data(DAT_FILE, DAT)