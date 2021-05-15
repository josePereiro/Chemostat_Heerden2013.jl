import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid, SpinLock

    #  ----------------------------------------------------------------------------
    # Run add https://github.com/josePereiro/Chemostat_Heerden2013.jl in the Julia Pkg REPL to install the
    # package, then you must activate the package enviroment (see README)
    import Chemostat_Heerden2013
    const ChHd = Chemostat_Heerden2013
    const Hd  = ChHd.HeerdenData;
    const BD  = ChHd.BegData;
    const iJR = ChHd.iJR904

    #  ----------------------------------------------------------------------------
    # run add "https://github.com/josePereiro/Chemostat" in the 
    # julia Pkg REPL for installing the package
    import Chemostat
    import Chemostat.LP.MathProgBase
    const Ch = Chemostat
    const ChU = Ch.Utils
    const ChLP = Ch.LP

    import JuMP, GLPK
    using Serialization
    import UtilsJL
    const UJL = UtilsJL
end

## -----------------------------------------------------------------------------------------------
const FBA_Z_FIX_MIN_COST        = :FBA_Z_FIX_MIN_COST
const FBA_Z_FIX_MAX_VG_MIN_COST = :FBA_Z_FIX_MAX_VG_MIN_COST
const FBA_MAX_Z_MIN_COST        = :FBA_MAX_Z_MIN_COST
const FBA_Z_FIX_MIN_VG_COST     = :FBA_Z_FIX_MIN_VG_COST
const FBA_Z_VG_FIX_MIN_COST     = :FBA_Z_VG_FIX_MIN_COST
const EXPS = 1:13

## -----------------------------------------------------------------------------------------------
# Data container
LPDAT = UJL.DictTree()

## -------------------------------------------------------------------
# FBA
let
    objider = iJR.BIOMASS_IDER
    costider = iJR.COST_IDER
    exglcider = iJR.EX_GLC_IDER
    max_sense = -1.0
    min_sense = 1.0

    iterator = Hd.val(:D) |> enumerate |> collect 
    for (exp, D) in iterator

        @info("Doing ", exp); println()

        # FBA_Z_FIX_MAX_VG_MIN_COST
        let
            model = iJR.load_model("fva_models", exp)
            exp_growth = Hd.val("D", exp)
            ChU.bounds!(model, objider, exp_growth, exp_growth)
            fbaout = ChLP.fba(model, exglcider, costider)

            LPDAT[FBA_Z_FIX_MAX_VG_MIN_COST, :model, exp] = model
            LPDAT[FBA_Z_FIX_MAX_VG_MIN_COST, :fbaout, exp] = fbaout
        end

        # FBA_MAX_Z_MIN_COST
        let
            model = iJR.load_model("fva_models", exp)
            fbaout = ChLP.fba(model, objider, costider)
            
            LPDAT[FBA_MAX_Z_MIN_COST, :model, exp] = model
            LPDAT[FBA_MAX_Z_MIN_COST, :fbaout, exp] = fbaout
        end

        # FBA_Z_FIX_MIN_COST
        let
            model = iJR.load_model("fva_models", exp)
            exp_growth = Hd.val("D", exp)
            ChU.bounds!(model, objider, exp_growth, exp_growth)
            fbaout = ChLP.fba(model, objider, costider)

            LPDAT[FBA_Z_FIX_MIN_COST, :model, exp] = model
            LPDAT[FBA_Z_FIX_MIN_COST, :fbaout, exp] = fbaout
        end

        # FBA_Z_FIX_MIN_VG_COST
        let
            model = iJR.load_model("fva_models", exp)
            exp_growth = Hd.val("D", exp)
            ChU.bounds!(model, objider, exp_growth, exp_growth)
            fbaout1 = ChLP.fba(model, exglcider; sense = max_sense)
            exglc = ChU.av(model, fbaout1, exglcider)
            ChU.bounds!(model, exglcider, exglc, exglc)
            fbaout = ChLP.fba(model, costider; sense = min_sense)

            LPDAT[FBA_Z_FIX_MIN_VG_COST, :model, exp] = model
            LPDAT[FBA_Z_FIX_MIN_VG_COST, :fbaout, exp] = fbaout
        end

        # FBA_Z_VG_FIX_MIN_COST
        let
            model = iJR.load_model("fva_models", exp)
            exp_growth = Hd.val("D", exp)
            ChU.bounds!(model, objider, exp_growth, exp_growth)
            exp_exglc = Hd.uval("GLC", exp)
            ChU.bounds!(model, exglcider, exp_exglc, exp_exglc)
            fbaout = ChLP.fba(model, costider; sense = min_sense)

            LPDAT[FBA_Z_VG_FIX_MIN_COST, :model, exp] = model
            LPDAT[FBA_Z_VG_FIX_MIN_COST, :fbaout, exp] = fbaout
        end
    end
end

## -------------------------------------------------------------------
LP_DAT_FILE = iJR.procdir("lp_dat_file.bson")
ChU.save_data(LP_DAT_FILE, LPDAT)

## -------------------------------------------------------------------
using Plots
let
    METHODS = [
        FBA_Z_FIX_MIN_COST, FBA_Z_FIX_MAX_VG_MIN_COST, 
        FBA_MAX_Z_MIN_COST, FBA_Z_FIX_MIN_VG_COST, 
        FBA_Z_VG_FIX_MIN_COST
    ]
    FLX_IDERS = ["GLC", "SUCC", "AC", "FORM"]
    rxn_map = iJR.load_rxns_map()
    color_pool = Plots.distinguishable_colors(length(FLX_IDERS))
    ider_color = Dict(ider => color for (ider, color) in zip(FLX_IDERS, color_pool))

    ps = Plots.Plot[]
    for method in METHODS
        p = plot(;title = string(method), xlabel = "exp", ylabel = "model")

        for exp in EXPS
            model = LPDAT[method, :model, exp]
            fbaout = LPDAT[method, :fbaout, exp]

            for Hd_ider in FLX_IDERS
                iJR_ider = rxn_map[Hd_ider]
                
                exp_val = Hd.uval(Hd_ider, exp) |> abs
                model_val = ChU.av(model, fbaout, iJR_ider) |> abs

                scatter!(p, [exp_val], [model_val]; label = "",
                    color = ider_color[Hd_ider], m = 8
                )
            end
        end
            push!(ps, p)
    end
    
    figname = string("fba_corrs")
    dir = iJR.plotsdir()
    UJL.mysavefig(ps, figname, dir)
end
