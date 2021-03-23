let
    method = ME_Z_EXPECTED_G_BOUNDED

    # Feed jobs
    Ch = Channel(1) do ch
        cGLCs = Hd.val("cGLC")
        for (exp, cGLC)  in enumerate(cGLCs)
            put!(ch, (exp, cGLC))
        end
    end

    @threads for thid in 1:nthreads()
        for (exp, cGLC) in Ch
            # handle cache
            datfile = dat_file(string(DAT_FILE_PREFFIX, method); exp)
            if isfile(datfile)
                lock(WLOCK) do
                    INDEX[method, :DFILE, exp] = datfile
                    @info("Cached loaded (skipping)",
                        exp, cGLC,datfile, thid
                    )
                    println()
                end
                continue
            end

            ## -------------------------------------------------------------------
            # setup
            model = base_model(exp)
            objidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
            M, N = size(model)
            exp_growth = Hd.val("D", exp)
            
            lock(WLOCK) do
                @info("Doing $(method)", exp, cGLC, thid)
                println()
            end
            
            ## -------------------------------------------------------------------
            # gradien descent
            epouts = Dict()
            x0 = 0.0
            @assert x0 == 0.0 # Must be zero for ME_Z_OPEN_G_OPEN
            x1 = 10.0
            maxΔ = 5e3
            th = 1e-3
            epout_seed = nothing
            target = exp_growth
            beta_vec = zeros(size(model, 2))

            upfrec_time = 10 # secunds
            last_uptime = time()
            gd_it = 1

            ## -------------------------------------------------------------------
            function upfun(beta)

                if haskey(epouts, beta) 
                    epout = epouts[beta]
                else
                    beta_vec[objidx] = beta
                    epout = ChEP.maxent_ep(model; 
                        beta_vec,
                        alpha = Inf,
                        damp = 0.9,
                        epsconv = 1e-4, 
                        verbose = false, 
                        solution = epout_seed,
                        maxiter = 5000
                    )
                    epouts[beta] = epout
                end
                epout_seed = epout
                ep_growth = ChU.av(epout)[objidx]

                update = gd_it == 1 || abs(last_uptime - time()) > upfrec_time || 
                    epout.status != ChEP.CONVERGED_STATUS
                update && lock(WLOCK) do
                    diff = abs.(exp_growth - ep_growth)
                    @info(
                        "Grad descent... ", 
                        exp, gd_it, 
                        epout.status, epout.iter, 
                        ep_growth, exp_growth, diff, 
                        beta,
                        thid
                    ); println()
                    gd_it += 1
                    last_uptime = time()
                end
                return ep_growth
            end

            ## -------------------------------------------------------------------
            # maxent
            expβ = UJL.grad_desc(upfun; x0, x1, th, maxΔ, 
                target, maxiters = 2000, verbose = false
            )

            ## -------------------------------------------------------------------
            lock(WLOCK) do
                # Storing
                dat = Dict()
                dat[:exp_beta] = expβ
                dat[:epouts] = epouts
                dat[:exp_beta] = maximum(keys(epouts))
                dat[:model] = model |> ChU.compressed_model

                # caching
                serialize(datfile, dat)
                INDEX[method, :DFILE, exp] = datfile

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