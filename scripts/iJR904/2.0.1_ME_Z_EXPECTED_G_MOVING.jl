let
    method = ME_Z_EXPECTED_G_MOVING

    # Feed jobs
    Ch = Channel(1) do ch
        cGLCs = Hd.val("cGLC")
        for (exp, cGLC)  in enumerate(cGLCs)
            put!(ch, (exp, cGLC))
        end
    end

    # @threads 
    for thid in 1:nthreads()
        for (exp, cGLC) in Ch

            ## -------------------------------------------------------------------
            # handle cache
            datfile = dat_file(;method, exp)
            check_cache(;method, exp) && continue

            ## -------------------------------------------------------------------
            # SetUp
            model_cache_id = (:MODEL0_CACHE, exp)
            # UJL.delete_cache(model_cache_id) # uncomment to reset
            model =  UJL.load_cache(model_cache_id; verbose = true) do
                BASE_MODELS = ChU.load_data(iJR.BASE_MODELS_FILE; verbose = false);
                model_dict = BASE_MODELS["fva_models"][exp]
                model0 = ChU.MetNet(;model_dict...) |> ChU.uncompressed_model
                ChU.scale_polytope!(model0, 0.5)
                return ChLP.fva_preprocess(model0; 
                    verbose = false, check_obj = iJR.BIOMASS_IDER
                )   
            end
            objidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
            M, N = size(model)
            exp_growth = Hd.val("D", exp)
            biom_lb, biom_ub = ChU.bounds(model, iJR.BIOMASS_IDER)
            if biom_ub < exp_growth
                lock(WLOCK) do
                    INDEX[method, :DFILE, exp] = :unfeasible
                    @info("Not feasible (skipping)", 
                        biom_ub ,exp_growth, 
                        thid
                    ); println()
                end
                continue
            end
            cgD_X = -Hd.cval(:GLC, exp) * Hd.val(:D, exp) / Hd.val(:DCW, exp)
            exglc_L = ChU.lb(model, iJR.EX_GLC_IDER)
            exglc_qta = abs(exglc_L * 0.005)
            expβ = 0.0
            epouts = Dict()
            epout = nothing

            for movround in 1:1000

                empty!(epouts)

                ## -------------------------------------------------------------------
                # GRAD DESCEND
                x0 = expβ
                x1 = 10.0
                maxΔ = max(expβ * 0.05, 5e3)
                gd_th = 1e-3
                target = exp_growth
                beta_vec = zeros(size(model, 2))
        
                upfrec_time = 10 # secunds
                last_uptime = time()
                gd_it = 1
        
                ## -------------------------------------------------------------------
                function upfun(beta)
        
                    beta_vec[objidx] = beta
                    epout = ChEP.maxent_ep(model; 
                        beta_vec,
                        alpha = Inf,
                        maxiter = 5000,
                        epsconv = 1e-3, 
                        verbose = false, 
                        solution = epout
                    )
                    epouts[beta] = epout
                    
                    # ep_growth = ChU.av(epout)[objidx]
                    ep_growth = ChU.av(model, epout, iJR.BIOMASS_IDER)
        
                    update = gd_it == 1 || abs(last_uptime - time()) > upfrec_time || 
                        epout.status != ChEP.CONVERGED_STATUS
        
                    update && lock(WLOCK) do
                        diff = abs.(exp_growth - ep_growth)
                        @info(
                            "Grad descent... ", 
                            exp, gd_it, 
                            epout.status, epout.iter, 
                            ep_growth, exp_growth, diff, 
                            beta, thid
                        ); println()
                        gd_it += 1
                        last_uptime = time()
                        
                    end
                    return ep_growth
                end
        
                ## -------------------------------------------------------------------
                # FIND BETA
                expβ = UJL.grad_desc(upfun; x0, x1, th = gd_th, maxΔ, 
                    target, maxiters = 2000, verbose = false
                )
        
                ## -------------------------------------------------------------------
                # MOVE V_UB
                Δstep = 0.5
                exglc_lb = ChU.lb(model, iJR.EX_GLC_IDER)
                exglc_ub = ChU.ub(model, iJR.EX_GLC_IDER)
                vg_avPME = ChU.av(model, epout, iJR.EX_GLC_IDER)

                # lb is the uptake limit
                dist = cgD_X - vg_avPME
                Δexglc_lb = sign(dist) * max(exglc_qta, abs(dist * Δstep))
                exglc_lb = min(exglc_ub,
                    max(exglc_L, exglc_lb + Δexglc_lb)
                )
                ChU.lb!(model, iJR.EX_GLC_IDER, exglc_lb)
                
                ## -------------------------------------------------------------------
                # INFO AND CONV
                ep_growth = ChU.av(model, epout, iJR.BIOMASS_IDER)
                gd_err = abs(exp_growth - ep_growth) / exp_growth
                conv = cgD_X <= vg_avPME && epout.status == ChEP.CONVERGED_STATUS && gd_err < gd_th
                
                lock(WLOCK) do
                    @info("Round Done", 
                        movround, conv,
                        gd_err, exglc_qta, Δexglc_lb,
                        (vg_avPME, cgD_X), 
                        exglc_ub, exglc_lb,  exglc_L, 
                        thid
                    ); println()
                end
                conv && break
        
            end #  for movround in 1:1000

            ## -------------------------------------------------------------------
            lock(WLOCK) do
                # Storing
                dat = Dict()
                dat[:exp_beta] = expβ
                dat[:epouts] = epouts
                dat[:model] = model |> ChU.compressed_model

                # caching
                serialize(datfile, dat)

                ep_growth = ChU.av(epouts[expβ])[objidx]
                diff = abs.(exp_growth - ep_growth)
                @info("Finished ",
                    exp, expβ, 
                    length(epouts),
                    ep_growth, exp_growth, diff, 
                    thid
                ); println()
            end

        end # for (exp, cGLC) in Ch
    end # for thid in 1:nthreads()
end

## -------------------------------------------------------------------
# Further convergence
let
    method = ME_Z_EXPECTED_G_MOVING
    objider = iJR.BIOMASS_IDER

    iterator = Hd.val("cGLC") |> enumerate |> collect 
    @threads for (exp, cGLC) in iterator

        datfile = dat_file(;method, exp)
        datfile == :unfeasible && continue
        dat = deserialize(datfile)
        model, epouts = ChU.uncompressed_model(dat[:model]) , dat[:epouts]

        exp_growth = Hd.val(:D, exp)
        exp_beta = dat[:exp_beta]
        exp_epout = epouts[exp_beta]

        lock(WLOCK) do
            @info("Converging...", 
                exp, method,
                exp_beta, exp_epout.status, 
                threadid()
            ); println()
        end
        converg_status = get!(dat, :converg_status, :undone)
        converg_status == :done && continue

        model = ChLP.fva_preprocess(model; verbose = false, 
            check_obj = iJR.BIOMASS_IDER
        )
        
        new_epout = nothing
        try
            objidx = ChU.rxnindex(model, objider)
            beta_vec = zeros(size(model, 2)); 
            beta_vec[objidx] = exp_beta
            new_epout = ChEP.maxent_ep(model; 
                beta_vec, alpha = Inf, 
                epsconv = 1e-5, verbose = false, 
                solution = exp_epout, maxiter = 5000
            )
        catch err; @warn("ERROR", err); println() end

        ep_growth = isnothing(new_epout) ? 0.0 : ChU.av(model, new_epout, objider)
        fail = isnan(ep_growth) || ep_growth == 0.0 
        epouts[exp_beta] = fail ? exp_epout : new_epout
        
        # Saving
        lock(WLOCK) do
            @info("Saving...", 
                exp, method, 
                exp_beta, 
                new_epout.status,
                new_epout.iter,
                threadid()
            ); println()
        end
        dat[:model] = ChU.compressed_model(model)
        dat[:converg_status] = :done
        serialize(datfile, dat)
    end
end