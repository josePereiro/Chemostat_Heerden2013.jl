let

    # Feed jobs
    Ch = Channel(1) do ch
        cGLCs = Hd.val("cGLC")
        for (exp, cGLC)  in enumerate(cGLCs)
            put!(ch, (exp, cGLC))
        end
    end
    
    @threads for thid in 1:nthreads()
        for (exp, cGLC) in Ch

            lock(WLOCK) do
                @info("Collecting $(ME_Z_OPEN_G_OPEN)", 
                    exp, cGLC, thid
                ); println()
            end

            exp_file = INDEX[ME_Z_EXPECTED_G_BOUNDED, :DFILE, exp]
            exp_dat = deserialize(exp_file)

            ME_Z_OPEN_G_OPEN_dat = Dict()
            ME_Z_OPEN_G_OPEN_dat[:exp_beta] = 0.0
            epout = exp_dat[:epouts][0.0]  # At beta 0
            ME_Z_OPEN_G_OPEN_dat[:epouts] = Dict(0.0 => epout)
            ME_Z_OPEN_G_OPEN_dat[:model] = exp_dat[:model]

            # save ME_Z_OPEN_G_OPEN
            ME_Z_OPEN_G_OPEN_file = dat_file(string(DAT_FILE_PREFFIX, ME_Z_OPEN_G_OPEN); exp)
            serialize(ME_Z_OPEN_G_OPEN_file, ME_Z_OPEN_G_OPEN_dat)
            INDEX[ME_Z_OPEN_G_OPEN, :DFILE, exp] = ME_Z_OPEN_G_OPEN_file
        end # for (exp, cGLC) in Ch
    end # @threads for thid in 1:nthreads()
end