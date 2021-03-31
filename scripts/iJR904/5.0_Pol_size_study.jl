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

    using Plots
    using ProgressMeter
    using Base.Threads
    using Serialization
    import UtilsJL
    const UJL = UtilsJL
end

##  ----------------------------------------------------------------------------
function load_model(modelkey::String = "max_model")
    BASE_MODELS = ChU.load_data(iJR.BASE_MODELS_FILE; verbose = false);
    model_dict = BASE_MODELS[modelkey]
    ChU.MetNet(;model_dict...) |> ChU.uncompressed_model
end

# -------------------------------------------------------------------
fileid = "5.0"
function mysavefig(p, pname; params...) 
    fname = UJL.mysavefig(p, string(fileid, "_", pname), iJR.MODEL_FIGURES_DIR; params...)
    @info "Plotting" fname
end

##  ----------------------------------------------------------------------------
let
    model0 = load_model("max_model")

    bins = 50
    Hd_rxns_map = iJR.load_Hd_rxns_map()
    offsetf = 1.1
    
    exglcidx = ChU.rxnindex(model0, iJR.EX_GLC_IDER)
    biomidx = ChU.rxnindex(model0, iJR.BIOMASS_IDER)
    exglcL, exglcU = ChU.bounds(model0, exglcidx)
    
    maxD = maximum(Hd.val(:D)) * offsetf 
    max_cgD_X = -maximum(Hd.ciD_X(:GLC)) * offsetf
    
    Ds = range(0.01, maxD; length = bins)
    cgD_Xs = range(max_cgD_X, exglcU; length = bins)

    box_vols = zeros(bins, bins) 
    
    Hd_ids = ["GLC", "SUCC", "AC", "FORM"]
    model_ids = [Hd_rxns_map[Hd_id] for Hd_id in Hd_ids]
    model_idxs = [ChU.rxnindex(model0, model_id) for model_id in model_ids]

    # feeding task
    Ch = Channel(nthreads()) do Ch_
        @showprogress for (Di, D) in enumerate(Ds)
            for (cgD_Xi, cgD_X) in enumerate(cgD_Xs)
                put!(Ch_, (Di, D, cgD_Xi, cgD_X))
            end
        end 
    end

    # compute volume map
    @threads for _ in 1:nthreads()
        thid = threadid()
        thid == 1 && continue
        model = deepcopy(model0)
        for (Di, D, cgD_Xi, cgD_X) in Ch
            
            # Reduce Pol
            ChU.lb!(model, exglcidx, cgD_X)
            ChU.bounds!(model, biomidx, D, D)
            try
                L, U = ChLP.fva(model, model_idxs)
                vol = prod(abs.(L .- U))
                box_vols[Di, cgD_Xi] = max(0.0, log10(vol + 1e-50))
            catch
                box_vols[Di, cgD_Xi] = NaN
            end
        end
    end

    # vol map
    p = heatmap(Ds, cgD_Xs, box_vols'; 
        title = "Polytope Box volume", label = "", 
        xlabel = "D", ylabel = "cgD_X"
    )

    # exp vals
    scatter!(p, Hd.val(:D), -Hd.ciD_X(:GLC);
        label = "exp", color = :white, 
        m = 8
    )

    mysavefig(p, "pol_box_volume"; bins)
end