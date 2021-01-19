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
end

## -------------------------------------------------------------------
# globals
const WLOCK = ReentrantLock()
#TODO make a proper caching system

## -------------------------------------------------------------------
# Prepare fva models
const DAT_FILE = iJR.MAXENT_FBA_EB_BOUNDLES_FILE
const DATA = isfile(DAT_FILE) ? ChU.load_data(DAT_FILE) : Dict()
let
    cGLCs = Hd.val("cGLC") |> enumerate |> collect
    for (exp, cGLC) in cGLCs 
        
        @info "Doing" exp cGLC length(DATA) threadid()
        println()

        D = get!(DATA, exp, Dict())
        haskey(D, "model0") && continue # check caching

        # prepare model
        model = ChU.MetNet(ChU.load_data(iJR.BASE_MODELS_FILE; verbose = true))
        model = ChU.uncompressed_model(model)
        objidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)

        exp_growth = Hd.val("D", exp)
        expξ = Hd.val("xi", exp)
        intake_info = deepcopy(iJR.base_intake_info)
        intake_info[iJR.EX_GLC_IDER]["c"] = cGLC
        ChSS.apply_bound!(model, expξ, intake_info)
        model = ChLP.fva_preprocess(model; check_obj = iJR.BIOMASS_IDER, 
            verbose = true)
        
        D["model0"] = ChU.compressed_model(model)
        
        # caching
        ChU.save_data(DAT_FILE, DATA; verbose = false)
    end
end

## -------------------------------------------------------------------
# Simulation
let
    
    cGLCs = Hd.val("cGLC") |> enumerate |> collect
    @threads for (exp, cGLC) in cGLCs 

        # handle cache
        D = DATA[exp]
        haskey(D, "epouts") && continue

        # setup
        model = D["model"]
        model = ChU.uncompressed_model(model)
        objidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
        M, N = size(model)
        exp_growth = Hd.val("D", exp)
        expξ = Hd.val("xi", exp)
        
        lock(WLOCK) do
            @info "Doing" exp cGLC threadid()
            println()
        end
        
        # gradien descent
        epouts = Dict()
        x0 = [0.0]
        x1 = [10.0]
        C = [5e3]
        th = 1e-3
        epout_seed = nothing
        target = [exp_growth]
        beta_vec = zeros(size(model, 2))

        upfrec = 10
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
                    epsconv = 1e-5, 
                    verbose = false, 
                    solution = epout_seed,
                    maxiter = 5000
                )
                epouts[beta] = epout
            end
            epout_seed = epout
            curr = [ChU.av(epout)[objidx]]

            if it == 1 || rem(it, upfrec) == 0
                lock(WLOCK) do
                    @info "Grad descent... " exp it beta abs.(target .- curr) threadid()
                    println()
                    it += 1
                end
            end
            curr
        end

        # maxent
        expβ = ChSU.grad_desc(upfun; x0, x1, th, C, 
            target, maxiters = 1500, verbose = false) |> first

        # fba
        fbaout = let
            lmodel = deepcopy(model)
            ChU.ub!(lmodel, iJR.BIOMASS_IDER, Hd.val("D", exp))
            ChLP.fba(lmodel, iJR.BIOMASS_IDER, iJR.COST_IDER)
        end

        lock(WLOCK) do
            # Storing
            D["exp_beta"] = expβ
            D["exp_xi"] = expξ
            D["epouts"] = epouts
            D["fbaout"] = fbaout
            D["epouts"] = epouts

            # caching
            ChU.save_data(DAT_FILE, DATA; verbose = false)

            @info "Finished " exp expβ length(epouts) threadid()
            println()
        end
    end
end

## -------------------------------------------------------------------
# Fva
# let
#     ϵ = 0.95 # relax factor
#     @info "Doing FVA" ϵ
#     for (exp, D) in DATA
        
#         model = D["model"]
#         fbaout = D["fbaout"]

#         biom_val = ChU.av(model, fbaout, iJR.BIOMASS_IDER)
#         cost_val = ChU.av(model, fbaout, iJR.COST_IDER)
#         biom_bounds = ChU.bounds(model, iJR.BIOMASS_IDER)
#         cost_bounds = ChU.bounds(model, iJR.COST_IDER)
        
#         @info "Doing" exp biom_bounds cost_bounds
#         ChU.bounds!(model, iJR.BIOMASS_IDER, ϵ * biom_val, biom_val)
#         ChU.bounds!(model, iJR.COST_IDER, ϵ * cost_val, cost_val)
        
#         D["fva"] = ChLP.fva(model)
        
#         ChU.bounds!(model, iJR.BIOMASS_IDER, biom_bounds...)
#         ChU.bounds!(model, iJR.COST_IDER, cost_bounds...)
#     end
# end

## -------------------------------------------------------------------
# save results
ChU.save_data(DAT_FILE, DATA; verbose = false)