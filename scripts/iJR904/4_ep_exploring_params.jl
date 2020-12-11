import DrWatson: quickactivate, savename
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

function partial_test(model, title  = "PARTIAL TEST")
    withcost = iJR.COST_IDER in model.rxns
    iders = withcost ? [iJR.BIOMASS_IDER, iJR.COST_IDER] : [iJR.BIOMASS_IDER]
    fbaout = ChLP.fba(model, iders...);
    ChU.tagprintln_inmw(title, 
        "\nsize:             ", size(model),
        "\nobj_ider:         " , iJR.BIOMASS_IDER,
        "\nfba obj_val:      ", ChU.av(model, fbaout, iJR.BIOMASS_IDER),
        "\ncost_ider:        ", withcost ? iJR.COST_IDER : "not found",
        "\nfba cost_val:     ", withcost ? ChU.av(model, fbaout, iJR.COST_IDER) : "not found",
    )
end

## ----------------------------------------------------------------------
function prepare_model(exp; verbose = true)
    # prepare model
    model = ChU.MetNet(ChU.load_data(iJR.BASE_MODEL_FILE; verbose))
    model = ChU.uncompressed_model(model)

    exp_growth = Hd.val("D", exp)
    expξ = Hd.val("xi", exp)
    intake_info = deepcopy(iJR.base_intake_info)
    intake_info["EX_glc_LPAREN_e_RPAREN_"]["c"] = Hd.val("cGLC", exp)
    ChSS.apply_bound!(model, expξ, intake_info)
    return model
end

## ----------------------------------------------------------------------
# # Data
let
    # setup
    exp = 1    
    # D = get!(DATA, "iJR904", Dict())
    model = prepare_model(exp)
    obj_ider = iJR.BIOMASS_IDER
    obj_idx = ChU.rxnindex(model, obj_ider)
    alpha = Inf
    epsconv = 1e-5
    maxiter = Int(1e7)
    damp = 0.99
    epouts = Dict()
    fbaout = ChLP.fba(model, obj_idx)
    partial_test(model); println()

    

    # seed
    seed_file = joinpath(iJR.MODEL_PROCESSED_DATA_DIR, "ept0_epout_seed.bson")
    if !isfile(seed_file)
        seed = ChEP.maxent_ep(model; alpha, epsconv, damp, maxiter)
        ChU.save_data(seed_file, seed)
    else
        seed = ChU.load_data(seed_file)
    end
    
    # Simulations
    write_lock = ReentrantLock()

    beta_range = [0.0; 10.0.^(-1:0.1:10)]
    var_orders = collect(5:25)
    
    @threads for var_order in var_orders

        lmodel = deepcopy(model)
        M, N = size(model)
        minvar, maxvar = 10.0^(-var_order), 10.0^(var_order)

        # MaxEnt EP
        beta_vec = zeros(N)
        epoutT0 = deepcopy(seed)

        D = Dict()
        D["epouts"] = Dict()
        for beta in beta_range
            beta_vec[obj_idx] = beta

            lock(write_lock) do
                @info "Doing" threadid() var_order beta
                println()
            end
            
            try
                epoutT0 = ChEP.maxent_ep(lmodel; alpha, beta_vec, 
                    epsconv, minvar, maxvar, solution = epoutT0,
                    verbose = false
                )
            catch err
                lock(write_lock) do
                    @info string("At threadid() = ", threadid())
                    @error err_str(err; max_len = 500)
                    println()
                end
                break
            end
            epoutT0.status != ChEP.CONVERGED_STATUS && break

            D["epouts"][beta] = epoutT0
        end

        # store and save
        D["beta_range"] = collect(keys(D["epouts"])) |> unique |> sort
        D["model"] = lmodel
        D["minvar"], D["maxvar"] = minvar, maxvar
        D["epsconv"] = epsconv
        D["alpha"] = alpha

        fname = savename("exploration_", (;var_order), "bson")
        fpath = joinpath(iJR.MODEL_PROCESSED_DATA_DIR, fname)
        ChU.save_data(fpath, D)

        lock(write_lock) do
            println()
        end

    end # @threads for var_order in var_orders
end
