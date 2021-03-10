import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

@time begin
    import SparseArrays
    using Base.Threads
    using Serialization

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

    import UtilsJL
    const UJL = UtilsJL
    UJL.set_cache_dir(iJR.CACHE_DIR)
    
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

# -------------------------------------------------------------------
function base_model(exp)
    BASE_MODELS = ChU.load_data(iJR.BASE_MODELS_FILE; verbose = false);
    model_dict = BASE_MODELS["fva_models"][exp]
    ChU.MetNet(;model_dict...) |> ChU.uncompressed_model
end

# -------------------------------------------------------------------
# METHOD VARIANTS
const ME_Z_OPEN_G_OPEN          = :ME_Z_OPEN_G_OPEN           # Do not use extra constraints
const ME_Z_EXPECTED_G_BOUNDED   = :ME_Z_EXPECTED_G_BOUNDED    # Match ME and Dy biom average and constraint av_ug
const ME_Z_EXPECTED_G_MOVING    = :ME_Z_EXPECTED_G_MOVING     # 
const ME_Z_FIXXED_G_BOUNDED     = :ME_Z_FIXXED_G_BOUNDED      # Fix biom around observed

## -------------------------------------------------------------------
# ME_Z_EXPECTED_G_MOVING
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
            cgD_X = -Hd.cval(:GLC, exp) * Hd.val(:D, exp) / Hd.val(:DCW, exp)
            exglc_L = ChU.lb(model, iJR.EX_GLC_IDER)
            exglc_qta = abs(exglc_L * 0.005)
            expβ = 0.0
            epouts_cid = (:EPOUTS_CACHE, exp)
            epouts = ChU.load_cache(epouts_cid, Dict(); verbose = false)
            epout_cid = (:EPOUT_CACHE, exp)
            epout = ChU.load_cache(epout_cid; verbose = false)
            for movround in 1:1000

                ## -------------------------------------------------------------------
                # GRAD DESCEND
                x0 = expβ
                x1 = 10.0
                maxΔ = 5e3
                th = 1e-3
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
        
                        try;
                            ChU.save_cache(epouts_cid, epouts; verbose = false)
                            ChU.save_cache(epout_cid, epout; verbose = false)
                        catch err
                            @warn("ERROR SAVING IGNORED", err)
                        end
                    end
                    return ep_growth
                end
        
                ## -------------------------------------------------------------------
                # FIND BETA
                expβ = UJL.grad_desc(upfun; x0, x1, th, maxΔ, 
                    target, maxiters = 2000, verbose = false
                )
        
                ## -------------------------------------------------------------------
                # MOVE V_UB
                Δstep = 0.5
                exglc_lb = ChU.lb(model, iJR.EX_GLC_IDER)
                exglc_ub = ChU.ub(model, iJR.EX_GLC_IDER)
                vg_avPME = ChU.av(model, epout, iJR.EX_GLC_IDER)
                # vg_avPME = exglc_lb * 0.8
                # lb is the uptake limit
                dist = cgD_X - vg_avPME
                Δexglc_lb = sign(dist) * max(exglc_qta, abs(dist * Δstep))
                exglc_lb = min(exglc_ub,
                    max(exglc_L, exglc_lb + Δexglc_lb)
                )
                ChU.lb!(model, iJR.EX_GLC_IDER, exglc_lb)
                
                ## -------------------------------------------------------------------
                # INFO AND CONV
                conv = cgD_X <= vg_avPME && epout.status == ChEP.CONVERGED_STATUS
                
                lock(WLOCK) do
                    @info("Round Done", 
                        movround, conv,
                        dist, exglc_qta, Δexglc_lb,
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

## -------------------------------------------------------------------
# ME_Z_EXPECTED_G_BOUNDED and ME_Z_OPEN_G_OPEN
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

## -------------------------------------------------------------------
# ME_Z_FIXXED_G_BOUNDED
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

## -------------------------------------------------------------------
# ME_Z_OPEN_G_OPEN
# It was computed in ME_Z_EXPECTED_G_BOUNDED
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

## -------------------------------------------------------------------
# SAVE INDEX
ChU.save_data(iJR.MAXENT_VARIANTS_INDEX_FILE, INDEX; verbose = false)