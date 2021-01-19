import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid

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
WLOCK = ReentrantLock()
SIM_GLOBAL_ID = "iJR904_MAXENT_VARIANTS"

# ChU.save_cache(res_id, (exp, dat); headline = "CATCHING RESULTS\n")
## -------------------------------------------------------------------
# Prepare fva models

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
# b0seed(exp) = ChU.load_data(iJR.MAXENT_B0SEEDS_FILE; verbose = false)[exp][:epoutb0]

## -------------------------------------------------------------------
const HOMO = :HOMO
const BOUNDED = :BOUNDED
const EXPECTED = :EXPECTED

## -------------------------------------------------------------------
# let
#     cGLCs = Hd.val("cGLC")[3:3] |> enumerate |> collect
#     @threads for (exp, cGLC) in cGLCs 
#         epout = b0seed(exp)
#         @info "seed " exp epout.status mean(epout.μ) mean(epout.σ) mean(epout.av) mean(epout.va)
#     end
# end
## -------------------------------------------------------------------
# @time let
#     exp = 3
#     model = base_model(exp)
#     damp = 0.985
#     epsconv = 1e-4
#     maxiter = 3000
#     alpha = Inf
#     solution = ChEP.maxent_ep(model; damp, alpha, epsconv, maxiter)
#     epout = ChEP.maxent_ep(model; damp, alpha, solution, epsconv, maxiter)
#     @info "seed " exp epout.status mean(epout.μ) mean(epout.σ) mean(epout.av) mean(epout.va)
# end;

## -------------------------------------------------------------------
# EXPECTED and HOMO
let
    cGLCs = Hd.val("cGLC") |> enumerate |> collect
    @threads for (exp, cGLC) in cGLCs 

        # handle cache
        dfile = dat_file(string("maxent_ep_boundle_", EXPECTED); exp)
        if isfile(dfile)
            lock(WLOCK) do
                INDEX[EXPECTED, :DFILE, exp] = dfile
                @info("Cached loaded (skipping)",
                    exp, cGLC,dfile, threadid()
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

        # fba
        fbaout = let
            lmodel = deepcopy(model)
            ChU.ub!(lmodel, iJR.BIOMASS_IDER, Hd.val("D", exp))
            ChLP.fba(lmodel, iJR.BIOMASS_IDER, iJR.COST_IDER)
        end

        lock(WLOCK) do
            # Storing
            dat = Dict()
            dat[:exp_beta] = expβ
            dat[:epouts] = 
            dat[:model] = model |> ChU.compressed_model
            dat[:fbaout] = fbaout

            # caching
            serialize(dfile, dat)
            INDEX[EXPECTED, :DFILE, exp] = dfile

            ep_growth = ChU.av(epouts[expβ])[objidx]
            diff = abs.(exp_growth - ep_growth)
            @info("Finished ",
                exp, expβ, 
                length(epouts),
                ep_growth, exp_growth, diff, 
                threadid()
            ); println()
        end
    end
end

## -------------------------------------------------------------------
# BOUNDED
let
    biomass_f = 0.01

    cGLCs = Hd.val("cGLC") |> enumerate |> collect
    @threads for (exp, cGLC) in cGLCs 

        # handle cache
        dfile = dat_file(string("maxent_ep_boundle_", BOUNDED); exp)
        if isfile(dfile)
            lock(WLOCK) do
                INDEX[BOUNDED, :DFILE, exp] = dfile
                @info("Cached loaded (skipping)",
                    exp, cGLC, dfile, threadid()
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
            dat[:epout] = epout
            dat[:model] = model |> ChU.compressed_model

            # caching
            serialize(dfile, dat)
            INDEX[BOUNDED, :DFILE, exp] = dfile

            @info("Finished ", exp, threadid())
            println()
        end
    end
end

## -------------------------------------------------------------------
# HOMO
# It was computed in EXPECTED
let
    cGLCs = Hd.val("cGLC") |> enumerate |> collect
    @threads for (exp, cGLC) in cGLCs 

        lock(WLOCK) do
            @info("Collecting $(HOMO)", 
                exp, cGLC, threadid()
            ); println()
        end

        exp_file = INDEX[EXPECTED, :DFILE, exp]
        exp_dat = deserialize(exp_file)

        homo_dat = Dict()
        homo_dat[:epout] = exp_dat[:epouts][0.0] # At beta 0
        homo_dat[:model] = exp_dat[:model]
        homo_dat[:fbaout] = exp_dat[:fbaout]

        # save homo
        homo_file = dat_file(string("maxent_ep_boundle_", HOMO); exp)
        serialize(homo_file, homo_dat)
        INDEX[HOMO, :DFILE, exp] = homo_file
            
    end
end

## -------------------------------------------------------------------
# SAVE INDEX
ChU.save_data(iJR.MAXENT_VARIANTS_INDEX_FILE, INDEX; verbose = false)