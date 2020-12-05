import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

import SparseArrays
import Base.Threads: @threads, threadid, SpinLock

## -------------------------------------------------------------------
# Run add https://github.com/josePereiro/Chemostat_Heerden2013.jl in the Julia Pkg REPL to install the
# package, then you must activate the package enviroment (see README)
import Chemostat_Heerden2013
const Hd  = Chemostat_Heerden2013.HeerdenData;
const BD  = Chemostat_Heerden2013.BegData;
const iJR = Chemostat_Heerden2013.iJR904

## -------------------------------------------------------------------
# run add "https://github.com/josePereiro/Chemostat" in the 
# julia Pkg REPL for installing the package
import Chemostat
import Chemostat.LP.MathProgBase
const ChU = Chemostat.Utils
const ChSS = Chemostat.SteadyState
const ChLP = Chemostat.LP
const ChEP = Chemostat.MaxEntEP
const ChSU = Chemostat.SimulationUtils
ChU.set_cache_dir(iJR.CACHE_DIR)

## -------------------------------------------------------------------
# Simulation
const DATA = Dict()
let
    verbose_lock = ReentrantLock()
    FREE = -1
    verbose_token = FREE
    
    cGLCs = Hd.val("cGLC")
    @threads for (exp, cGLC) in cGLCs |> enumerate |> collect

        # handle cache
        haskey(DATA, exp) && continue

        # handle verbose
        imverbose = false
        lock(verbose_lock) do
            if verbose_token == FREE
                verbose_token = threadid()
                imverbose = true
            end
        end

        # model
        model = ChU.MetNet(ChU.load_data(iJR.BASE_MODEL_FILE; verbose = imverbose))
        model = ChU.uncompressed_model(model)
        objidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
        M, N = size(model)
            
        imverbose && @info "Doing" threadid() exp cGLC 
        D = get!(DATA, exp, Dict())

        # prepare model
        exp_growth = Hd.val("D", exp)
        expξ = Hd.val("xi", exp)
        intake_info = deepcopy(iJR.base_intake_info)
        intake_info["EX_glc_LPAREN_e_RPAREN_"]["c"] = cGLC
        ChSS.apply_bound!(model, expξ, intake_info)
        
        # gradien descent
        epouts = Dict()
        x0 = [0.0]
        x1 = [10.0]
        C = [5e3]
        th = 1e-4
        epout_seed = nothing
        target = [Hd.val("D", exp)]
        beta_vec = zeros(size(model, 2))

        function upfun(betas)
            beta = first(betas)
            if haskey(epouts, beta) 
                epout = epouts[beta]
            else
                beta_vec[objidx] = beta
                epout = ChEP.maxent_ep(model; 
                    beta_vec,
                    alpha = Inf,
                    damp = 0.98,
                    epsconv = 1e-4, 
                    verbose = imverbose && isempty(epouts), 
                    solution = epout_seed,
                    maxiter = 3000
                )
                epouts[beta] = epout
            end
            epout_seed = epout
            r = [ChU.av(epout)[objidx]]
        end

        expβ = ChSU.grad_desc(upfun; x0, x1, th, C, 
            target, maxiters = 1500, verbose = imverbose)

        D["model"] = ChU.compressed_model(model)
        D["exp_beta"] = first(expβ)
        D["exp_xi"] = expξ
        D["epouts"] = epouts
        D["fbaout"] = ChLP.fba(model, iJR.BIOMASS_IDER, iJR.COST_IDER)
        D["epouts"] = epouts
        
        # release token
        lock(verbose_lock) do
            imverbose && (verbose_token = FREE; println())
        end
    end
end

## -------------------------------------------------------------------
# save results
ChU.save_data(iJR.MAXENT_FBA_EB_BOUNDLES_FILE, DATA)