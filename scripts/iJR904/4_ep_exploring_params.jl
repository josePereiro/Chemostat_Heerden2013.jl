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

## ----------------------------------------------------------------------
# Tooling
function err_str(err; max_len = 500)
    s = sprint(showerror, err, catch_backtrace())
    return length(s) > max_len ? s[1:max_len] * "\n[...]" : s
end

function prepare_model()
    model = ChU.MetNet(ChU.load_data(iJR.BASE_MODEL_FILE)) 

end

## -------------------------------------------------------------------
# Tools
function partial_test(model, title  = "PARTIAL TEST")
    withcost = iJR.COST_IDER in model.rxns
    iders = withcost ? [iJR.BIOMASS_IDER, iJR.COST_IDER] : [iJR.BIOMASS_IDER]
    fbaout = ChLP.fba(model, iders...);
    ChU.tagprintln_inmw(title, 
        "\nsize:             ", size(model),
        "\nobj_ider:         " , iJR.BIOMASS_IDER,
        "\nfba obj_val:      ", ChU.av(model, fbaout, iJR.BIOMASS_IDER),
        "\nmax exp obj_val:  ", maximum(Hd.val("D")),
        "\ncost_ider:        ", withcost ? iJR.COST_IDER : "not found",
        "\nfba cost_val:     ", withcost ? ChU.av(model, fbaout, iJR.COST_IDER) : "not found",
    )
end

my_try(f; max_len = 500) = try; (f(); return true) catch err; 
    (@error(err_str(err; max_len)); return false) end

## ----------------------------------------------------------------------
# Data
DATA = Dict()
let

    D = get!(DATA, "EColi_Core_FVA", Dict())
    D["model"] = ChT.ecoli_core_model();
    ChU.bounds!(D["model"], "SUCDi", 0.0, 0.0)
    D["model"] = ChLP.fva_preprocess(D["model"])
    D["obj_ider"] = ChT.ECOLI_MODEL_BIOMASS_IDER
    D["obj_idx"] = ChU.rxnindex(D["model"], D["obj_ider"])
    D["T!0alpha"] = 1e7
    D["T0alpha"] = Inf
    D["beta_range"] = [0.0; 10.0.^(-1:0.1:10)]
    
    # Simulations
    for (modelid, D) in DATA
        
        model = D["model"]
        M, N = size(model)

        # FBA
        D["fbaout"] = ChLP.fba(model, D["obj_idx"])

        # MaxEnt EP
        D["epouts"] = Dict()
        minvar, maxvar = 1e-18, 1e18
        maxiter = Int(1e6)
        verbose = true
        epsconv = 1e-5
        beta_vec = zeros(N)
        epoutT!0, epoutT0 = nothing, nothing

        for beta in D["beta_range"]
            beta_vec[D["obj_idx"]] = beta
            @info "" modelid beta
            
            # T>0
            my_try() do
                epoutT!0 = ChEP.maxent_ep(model; alpha = D["T!0alpha"], beta_vec, epsconv,
                    minvar, maxvar, maxiter, solution = epoutT!0, verbose);
            end && epoutT!0.status == ChEP.CONVERGED_STATUS || break

            # T0
            my_try() do
                epoutT0 = ChEP.maxent_ep(model; alpha = D["T0alpha"], beta_vec, epsconv,
                    minvar, maxvar, maxiter, solution = epoutT0, verbose);
            end && epoutT0.status == ChEP.CONVERGED_STATUS || break

            D["epouts"][(D["T!0alpha"], beta)] = epoutT!0
            D["epouts"][(D["T0alpha"], beta)] = epoutT0
        end

        D["beta_range"] = last.(keys(D["epouts"])) |> unique |> sort
    end
end;
