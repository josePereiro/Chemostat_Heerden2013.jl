import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

## -------------------------------------------------------------------
# Run add https://github.com/josePereiro/Chemostat_Heerden2013.jl in the Julia Pkg REPL to install the
# package, then you must activate the package enviroment (see README)
import Chemostat_Heerden2013
const Hd  = Chemostat_Heerden2013.HeerdenData;
const ECC = Chemostat_Heerden2013.EColiCore

## -------------------------------------------------------------------
# run add "https://github.com/josePereiro/Chemostat" in the 
# julia Pkg REPL for installing the package
import Chemostat
const ChU = Chemostat.Utils
const ChSS = Chemostat.SteadyState
const ChLP = Chemostat.LP
const ChEP = Chemostat.MaxEntEP
const ChSU = Chemostat.SimulationUtils

## -------------------------------------------------------------------
# Simulation
DATA = Dict()
let
    model = ChU.MetNet(ChU.load_data(ECC.BASE_MODEL_FILE))
    objidx = ChU.rxnindex(model, ECC.BIOMASS_IDER)
    M, N = size(model)

    for (exp, cGLC) in enumerate(Hd.val("cGLC"))
        
        @info "Doing" exp cGLC
        D = get!(DATA, exp, Dict())

        # prepare model
        exp_growth = Hd.val("D", exp)
        expξ = Hd.val("xi", exp)
        intake_info = deepcopy(ECC.base_intake_info)
        intake_info["EX_glc__D_e"]["c"] = cGLC
        ChSS.apply_bound!(model, expξ, intake_info)
        # ChU.bounds!(model, ECC.BIOMASS_IDER, 0.0, exp_growth * 1.1)

        # maxent ep params
        alpha = Inf
        beta_vec = zeros(size(model, 2))
        epsconv = 1e-7
        alpha = Inf
        verbose = false

        # gradien descent
        epouts = Dict()
        x0 = [0.0]
        x1 = [10.0]
        C = [5e2]
        th = 1e-4
        epout_seed = nothing
        target = [Hd.val("D", exp)]
        expβ = ChSU.grad_desc(;x0, x1, th, C, target, maxiters = 1500) do betas
            beta = first(betas)
            if haskey(epouts, beta) 
                epout = epouts[beta]
            else
                beta_vec[objidx] = beta
                solution = epout_seed
                epout = ChEP.maxent_ep(model; alpha, beta_vec, epsconv, verbose, solution)
                epouts[beta] = epout
            end
            epout_seed = epout
            r = [ChU.av(epout)[objidx]]
        end

        # collect results
        D["model"] = model
        D["exp_beta"] = first(expβ)
        D["exp_xi"] = expξ
        D["epouts"] = epouts
        D["fbaout"] = ChLP.fba(model, ECC.BIOMASS_IDER)

        println()
    end
end

## -------------------------------------------------------------------
# save results
ChU.save_data(ECC.MAXENT_RES_FILE, DATA)