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

# ----------------------------------------------------------------------------------
FLX_IDERS = ["GLC", "SUCC", "AC", "FORM"]
DAT[:FLX_IDERS] = FLX_IDERS;

# -------------------------------------------------------------------
# ME methods
const ME_Z_OPEN_G_OPEN        = :ME_Z_OPEN_G_OPEN
const ME_MAX_POL              = :ME_MAX_POL
const ME_MAX_POL_B0           = :ME_MAX_POL_B0
const ME_Z_EXPECTED_G_MOVING  = :ME_Z_EXPECTED_G_MOVING
const ME_Z_EXPECTED_G_BOUNDED = :ME_Z_EXPECTED_G_BOUNDED
const ME_Z_FIXXED_G_BOUNDED   = :ME_Z_FIXXED_G_BOUNDED

# LP methods
const FBA_Z_FIX_MIN_COST      = :FBA_Z_FIX_MIN_COST
const FBA_MAX_Z_MIN_COST   = :FBA_MAX_Z_MIN_COST
const FBA_Z_FIX_MIN_VG_COST   = :FBA_Z_FIX_MIN_VG_COST
const FBA_Z_VG_FIX_MIN_COST   = :FBA_Z_VG_FIX_MIN_COST

LP_METHODS = [
    FBA_Z_FIX_MIN_COST, FBA_MAX_Z_MIN_COST, 
    FBA_Z_FIX_MIN_VG_COST, FBA_Z_VG_FIX_MIN_COST
]
DAT[:LP_METHODS] = LP_METHODS

ME_METHODS = [
    # ME_Z_OPEN_G_OPEN, 
    ME_MAX_POL,
    ME_MAX_POL_B0,
    # ME_Z_FIXXED_G_BOUNDED, 
    ME_Z_EXPECTED_G_BOUNDED, 
    # ME_Z_EXPECTED_G_MOVING
]
DAT[:ME_METHODS] = ME_METHODS

ALL_METHODS = [LP_METHODS; ME_METHODS]
DAT[:ALL_METHODS] = ALL_METHODS

EXPS = Hd.EXPS 
DAT[:EXPS] = EXPS

# -------------------------------------------------------------------
ME_INDEX_FILE = iJR.procdir("maxent_ep_index.bson")
ME_INDEX = ChU.load_data(ME_INDEX_FILE; verbose = false);

# -------------------------------------------------------------------
LP_DAT_FILE = iJR.procdir("lp_dat_file.bson")
LP_DAT = ChU.load_data(LP_DAT_FILE; verbose = false);

# ----------------------------------------------------------------------------------
Hd_mets_map = iJR.load_mets_map()
Hd_rxns_map = iJR.load_rxns_map()

# ----------------------------------------------------------------------------------
const DAT_FILE_PREFFIX =  "maxent_ep_dat"

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

    # util fun
    isexpdep(method) = (method != ME_MAX_POL_B0)
    depks(method, typ, Kd_met, exp) = 
        isexpdep(method) ? (method, typ, Kd_met, exp) : (method, typ, Kd_met)
    function dat_file(;method, exp)
        kwargs = isexpdep(method) ? (;method, exp) : (;method)
        fname = UJL.mysavename(DAT_FILE_PREFFIX, "jls"; kwargs...)
        iJR.procdir(fname)
    end

    # Feed jobs
    nths = nthreads()
    Ch = Channel(nths) do ch
        for method in DAT[:ME_METHODS]
            for exp in DAT[:EXPS]
                put!(ch, (exp, method))
                !isexpdep(method) && break # Do only once
            end
        end
    end

    @threads for thid in 1:nths
        thid = threadid()
        thid == 1 && nths > 1 && continue
        for (exp, method) in Ch
            
            # ME data
            datfile = dat_file(;method, exp)
            !isfile(datfile) && continue
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
            
            lock(WLOCK) do
                DAT[depks(method, :flx, "D", exp)...] = ep_biom
                DAT[depks(method, :err, "D", exp)...] = ep_std
            end

            # exchanges
            for Hd_met in FLX_IDERS
                model_exch = Hd_rxns_map[Hd_met]

                # flxs
                ep_av = ChU.av(model, epout, model_exch)
                ep_std = sqrt(ChU.va(model, epout, model_exch))

                # proj 2d
                proj = ChLP.projection2D(model, objider, model_exch; l = 50)
                        
                lock(WLOCK) do
                    DAT[depks(method, :proj, Hd_met, exp)...] = proj
                    DAT[depks(method, :flx, Hd_met, exp)...] = ep_av
                    DAT[depks(method, :err, Hd_met, exp)...] = ep_std
                end
            end

             # additional fluxs
             for (ider, model_iders) in iJR.load_inners_idermap()
                # flxs
                ep_av = ChU.av(model, epout, model_iders[1])
                ep_std = sqrt(ChU.va(model, epout, model_iders[1]))
                if length(model_iders) == 2 # reversible
                    # r = r+ - r-
                    ep_av -= ChU.av(model, epout, model_iders[2])
                    ep_std += sqrt(ChU.va(model, epout, model_iders[2]))
                end

                # proj 2d (fwd only)
                proj = ChLP.projection2D(model, objider, model_iders[1]; l = 50)
                        
                lock(WLOCK) do
                    DAT[depks(method, :proj, ider, exp)...] = proj
                    DAT[depks(method, :flx, ider, exp)...] = ep_av
                    DAT[depks(method, :err, ider, exp)...] = ep_std
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

            # exchanges
            for Hd_ider in FLX_IDERS
                model_ider = Hd_rxns_map[Hd_ider]

                Hd_flx = Hd.uval(Hd_ider, exp)
                fba_flx = ChU.av(model, fbaout, model_ider)
                DAT[method, :flx, Hd_ider, exp] = fba_flx
            end

            # additional fluxs
            for (ider, model_iders) in iJR.load_inners_idermap()
                # flxs
                fba_flx = ChU.av(model, fbaout, model_iders[1])
                if length(model_iders) == 2 # reversible
                    # r = r+ - r-
                    fba_flx -= ChU.av(model, fbaout, model_iders[2])
                end
                        
                DAT[method, :flx, ider, exp] = fba_flx
            end

        end

    end # for method
end

## ----------------------------------------------------------------------------------
DAT_FILE = iJR.procdir("dat.bson")
UJL.save_data(DAT_FILE, DAT)