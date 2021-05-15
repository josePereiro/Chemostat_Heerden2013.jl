
##  ----------------------------------------------------------------------------
# pol vox volumen
let
    model0 = iJR.load_model("max_model")

    bins = 50
    Hd_rxns_map = iJR.load_rxns_map()
    offsetf = 1.1
    
    exglcidx = ChU.rxnindex(model0, iJR.EX_GLC_IDER)
    biomidx = ChU.rxnindex(model0, iJR.BIOMASS_IDER)
    exglcL, exglcU = ChU.bounds(model0, exglcidx)
    
    maxD = maximum(Hd.val(:D)) * offsetf 
    max_cgD_X = -maximum(Hd.ciD_X(:GLC)) * offsetf
    
    Ds = range(0.01, maxD; length = bins)
    cgD_Xs = range(max_cgD_X, exglcU; length = bins)

    box_vols = zeros(bins, bins) 
    
    Hd_ids = ["GLC", "AC", "NH4"]
    model_ids = [Hd_rxns_map[Hd_id] for Hd_id in Hd_ids]
    model_idxs = [ChU.rxnindex(model0, model_id) for model_id in model_ids]

    # feeding task
    nths = nthreads()
    Ch = Channel(nths) do Ch_
        @showprogress for (Di, D) in enumerate(Ds)
            for (cgD_Xi, cgD_X) in enumerate(cgD_Xs)
                put!(Ch_, (Di, D, cgD_Xi, cgD_X))
            end
        end 
    end

    # compute volume map
    @threads for _ in 1:nths
        thid = threadid()
        thid == 1 && nths > 1 && continue
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
        xlabel = "D", ylabel = "cgD/ X"
    )

    # exp vals
    exp_Ds = [Hd.val(:D, exp) for exp in EXPS]
    exp_cgD_Xs = [-Hd.ciD_X(:GLC, exp) for exp in EXPS]
    scatter!(p, exp_Ds, exp_cgD_Xs;
        label = "exp", color = :white, 
        m = 8
    )

    mysavefig(p, "pol_box_volume"; bins)
end