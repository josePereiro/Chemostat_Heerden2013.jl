import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

@time begin
    import SparseArrays
    using Base.Threads

    # -------------------------------------------------------------------
    # Run add https://github.com/josePereiro/Chemostat_Heerden2013.jl in 
    # the Julia Pkg REPL to install the package, then you must activate 
    # the package enviroment (see README)
    import Chemostat_Heerden2013
    const ChHd = Chemostat_Heerden2013
    const Hd  = ChHd.HeerdenData;
    const BD  = ChHd.BegData;
    const iJR = ChHd.iJR904

    # -------------------------------------------------------------------
    # run add "https://github.com/josePereiro/Chemostat" in the 
    # julia Pkg REPL for installing the package
    import Chemostat
    import Chemostat.LP.MathProgBase
    const Ch = Chemostat
    const ChU = Ch.Utils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP
    const ChEP = Ch.MaxEntEP
    const ChSU = Ch.SimulationUtils

    using Serialization
    # using Statistics
end

## -------------------------------------------------------------------
# globals
const WLOCK = ReentrantLock()
const SIM_GLOBAL_ID = "iJR904_MAXENT_VARIANTS"
const DAT_FILE_PREFFIX =  "maxent_ep_boundle_"

const INDEX = ChU.DictTree()
function dat_file(name; kwargs...)
    fname = ChU.mysavename(name, "jls"; kwargs...)
    joinpath(iJR.MODEL_PROCESSED_DATA_DIR, fname)
end

## -------------------------------------------------------------------
function base_model(exp)
    BASE_MODELS = ChU.load_data(iJR.BASE_MODELS_FILE; verbose = false);
    model_dict = BASE_MODELS["fva_models"][exp]
    ChU.MetNet(;model_dict...) |> ChU.uncompressed_model
end

## -------------------------------------------------------------------
const HOMO = :HOMO
const BOUNDED = :BOUNDED
const EXPECTED = :EXPECTED

## -------------------------------------------------------------------
# EXPECTED and HOMO
let

    Ch = Channel(1) do ch
        cGLCs = Hd.val("cGLC")
        for (exp, cGLC)  in enumerate(cGLCs)
            put!(ch, (exp, cGLC))
        end
    end

    @threads for thid in 1:nthreads()
        for (exp, cGLC) in Ch
            # handle cache
            datfile = dat_file(string(DAT_FILE_PREFFIX, EXPECTED); exp)
            if isfile(datfile)
                lock(WLOCK) do
                    INDEX[EXPECTED, :DFILE, exp] = datfile
                    @info("Cached loaded (skipping)",
                        exp, cGLC,datfile, threadid()
                    )
                    println()
                end
                continue
            end

            # setup
            model = base_model(exp)
            objidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
            M, N = size(model)
            exp_growth = Hd.val("D", exp)
            
            lock(WLOCK) do
                @info("Doing $(EXPECTED)", exp, cGLC, threadid())
                println()
            end
            
            # gradien descent
            epouts = Dict()
            x0 = [0.0] 
            @assert x0 == [0.0] # Must be zero for HOMO
            x1 = [10.0]
            C = [5e3]
            th = 1e-3
            epout_seed = nothing
            target = [exp_growth]
            beta_vec = zeros(size(model, 2))

            upfrec_time = 10 # secunds
            last_uptime = time()
            it = 1

            function upfun(betas)
                beta = first(betas)

                if haskey(epouts, beta) 
                    epout = epouts[beta]
                else
                    beta_vec[objidx] = beta
                    epout = ChEP.maxent_ep(model; 
                        beta_vec,
                        alpha = Inf,
                        damp = 0.985,
                        epsconv = 1e-4, 
                        verbose = false, 
                        solution = epout_seed,
                        maxiter = 5000
                    )
                    epouts[beta] = epout
                end
                epout_seed = epout
                ep_growth = ChU.av(epout)[objidx]

                update = it == 1 || abs(last_uptime - time()) > upfrec_time || 
                    epout.status != ChEP.CONVERGED_STATUS
                if update
                    lock(WLOCK) do
                        diff = abs.(exp_growth - ep_growth)
                        @info(
                            "Grad descent... ", 
                            exp, it, 
                            epout.status, epout.iter, 
                            ep_growth, exp_growth, diff, 
                            beta,
                            threadid()
                        ); println()
                        it += 1
                        last_uptime = time()
                    end
                end
                [ep_growth]
            end

            # maxent
            expβ = ChSU.grad_desc(upfun; x0, x1, th, C, 
                target, maxiters = 2000, verbose = false) |> first

            lock(WLOCK) do
                # Storing
                dat = Dict()
                dat[:exp_beta] = expβ
                dat[:epouts] = epouts
                dat[:model] = model |> ChU.compressed_model

                # caching
                serialize(datfile, dat)
                INDEX[EXPECTED, :DFILE, exp] = datfile

                ep_growth = ChU.av(epouts[expβ])[objidx]
                diff = abs.(exp_growth - ep_growth)
                @info("Finished ",
                    exp, expβ, 
                    length(epouts),
                    ep_growth, exp_growth, diff, 
                    threadid()
                ); println()
            end
        end # for (exp, cGLC) in Ch
    end # for thid in 1:nthreads()
end

## -------------------------------------------------------------------
# BOUNDED
let
    biomass_f = 0.01

    Ch = Channel(1) do ch
        cGLCs = Hd.val("cGLC")
        for (exp, cGLC)  in enumerate(cGLCs)
            put!(ch, (exp, cGLC))
        end
    end
    
    @threads for thid in 1:nthreads()
        for (exp, cGLC) in Ch
            # handle cache
            datfile = dat_file(string(DAT_FILE_PREFFIX, BOUNDED); exp)
            if isfile(datfile)
                lock(WLOCK) do
                    INDEX[BOUNDED, :DFILE, exp] = datfile
                    @info("Cached loaded (skipping)",
                        exp, cGLC, datfile, threadid()
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
                @info("Doing $(BOUNDED)", 
                    exp, cGLC, threadid()
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
                INDEX[BOUNDED, :DFILE, exp] = datfile

                @info("Finished ", exp, threadid())
                println()
            end
        end # for (exp, cGLC) in Ch
    end # for thid in 1:nthreads()
end

## -------------------------------------------------------------------
# HOMO
# It was computed in EXPECTED
let
    Ch = Channel(1) do ch
        cGLCs = Hd.val("cGLC")
        for (exp, cGLC)  in enumerate(cGLCs)
            put!(ch, (exp, cGLC))
        end
    end
    
    @threads for thid in 1:nthreads()
        for (exp, cGLC) in Ch

            lock(WLOCK) do
                @info("Collecting $(HOMO)", 
                    exp, cGLC, threadid()
                ); println()
            end

            exp_file = INDEX[EXPECTED, :DFILE, exp]
            exp_dat = deserialize(exp_file)

            homo_dat = Dict()
            homo_dat[:exp_beta] = 0.0
            epout = exp_dat[:epouts][0.0]  # At beta 0
            homo_dat[:epouts] = Dict(0.0 => epout)
            homo_dat[:model] = exp_dat[:model]

            # save homo
            homo_file = dat_file(string(DAT_FILE_PREFFIX, HOMO); exp)
            serialize(homo_file, homo_dat)
            INDEX[HOMO, :DFILE, exp] = homo_file
        end
    end
end

## -------------------------------------------------------------------
# SAVE INDEX
ChU.save_data(iJR.MAXENT_VARIANTS_INDEX_FILE, INDEX; verbose = false)