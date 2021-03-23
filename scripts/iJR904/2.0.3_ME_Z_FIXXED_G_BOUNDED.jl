let
    biomass_f = 0.01

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
            datfile = dat_file(string(DAT_FILE_PREFFIX, ME_Z_FIXXED_G_BOUNDED); exp)
            if isfile(datfile)
                lock(WLOCK) do
                    INDEX[ME_Z_FIXXED_G_BOUNDED, :DFILE, exp] = datfile
                    @info("Cached loaded (skipping)",
                        exp, cGLC, datfile, thid
                    ); println()
                end
                continue
            end

            # setup
            model = base_model(exp)
            objidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
            M, N = size(model)
            exp_growth = Hd.val("D", exp)
            ChU.ub!(model, iJR.BIOMASS_IDER, exp_growth * (1.0 + biomass_f))
            ChU.lb!(model, iJR.BIOMASS_IDER, exp_growth * (1.0 - biomass_f))
            model = ChLP.fva_preprocess(model, 
                batchlen = 50,
                check_obj = iJR.BIOMASS_IDER,
                verbose = true
            )

            lock(WLOCK) do
                @info("Doing $(ME_Z_FIXXED_G_BOUNDED)", 
                    exp, cGLC, thid
                ); println()
            end

            # maxent
            epout = ChEP.maxent_ep(model; alpha = Inf, damp = 0.985, epsconv = 1e-4, 
                        verbose = false, maxiter = 5000
                    )
                
            # storing
            lock(WLOCK) do
                # Storing
                dat = Dict()
                dat[:exp_beta] = 0.0
                dat[:epouts] = Dict(0.0 => epout)
                dat[:model] = model |> ChU.compressed_model

                # caching
                serialize(datfile, dat)
                INDEX[ME_Z_FIXXED_G_BOUNDED, :DFILE, exp] = datfile

                @info("Finished ", exp, thid)
                println()
            end
        end # for (exp, cGLC) in Ch
    end # for thid in 1:nthreads()
end