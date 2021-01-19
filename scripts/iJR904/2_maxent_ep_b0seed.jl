import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

# ------------------------------------------------------------------
# ARGS
using ArgParse

set = ArgParseSettings()
@add_arg_table! set begin
    "-w", "--workers"
        help = "number of workers to use"
        default = "1"
    "--init-clear"
        help = "clear cache before running the simulation"   
        action = :store_true
    "--finish-clear"
        help = "clear cache at the end"   
        action = :store_true
end

if isinteractive()
    # Dev vals
    wcount = 1
    init_clear_flag = false
    finish_clear_flag = false
else
    parsed_args = parse_args(set)
    wcount = parse(Int, parsed_args["workers"])
    init_clear_flag = parsed_args["init-clear"]
    finish_clear_flag = parsed_args["finish-clear"]
end

using Distributed
# -------------------------------------------------------------------
const NO_WORKERS = wcount
length(workers()) < NO_WORKERS && 
addprocs(NO_WORKERS; exeflags = "--project")
println("Working in: ", workers())


# -------------------------------------------------------------------
@time begin
    # Precompile in master
    import Chemostat_Heerden2013
    import Chemostat

    @everywhere begin
        
        import DrWatson: quickactivate 
        quickactivate(@__DIR__, "Chemostat_Heerden2013")

        using SparseArrays
        
        ## ----------------------------------------------------------------------------
        import Chemostat
        const Ch = Chemostat
        const ChU = Ch.Utils
        const ChSS = Ch.SteadyState
        const ChLP = Ch.LP
        const ChEP = Ch.MaxEntEP
        const ChSU = Ch.SimulationUtils

        ## ----------------------------------------------------------------------------
        import Chemostat_Heerden2013
        const ChH = Chemostat_Heerden2013
        const Hd  = ChH.HeerdenData;
        const BD  = ChH.BegData;
        const iJR = ChH.iJR904
        ChU.set_cache_dir(iJR.CACHE_DIR)

    end
end

## ------------------------------------------------------------------
# CLEAR CACHE (WARNING)
if init_clear_flag
    ChU.tagprintln_inmw("CLEARING CACHE ")
    ChU.delete_temp_caches()
    ChU.println_inmw("\n")
end

## ------------------------------------------------------------------
# SIMULATION GLOBAL ID
# This must uniquely identify this simulation version
# It is used to avoid cache collisions
@everywhere const SIM_GLOBAL_ID = "iJR904_MAXENT_B0SEED"

## ----------------------------------------------------------------------------
# SIM IDS
# Collect all the computed simulation ids for collecting results
const REMCHNL = RemoteChannel(() -> Channel{Any}(10))
const RES_IDS = []
const RESCOLL = @async while true
    id = take!(REMCHNL)
    push!(RES_IDS, id)
end

## ----------------------------------------------------------------------------
@everywhere function prepare_model(exp)
    BASE_MODELS = ChU.load_data(iJR.BASE_MODELS_FILE; verbose = false);
    model_dict = BASE_MODELS["fva_models"][exp]
    ChU.MetNet(;model_dict...) 
end

## ----------------------------------------------------------------------------
# RUN SIMULATION
empty!(RES_IDS)
cGLCs = Hd.val("cGLC") # Test
pmap(enumerate(cGLCs)) do (exp, cGLC) # This is parallizable

    EPOCHLEN = 100
    ξ = Hd.val("xi", exp)
    β0 = 0.0

    ChU.tagprintln_inmw("PROCESSING EXPERIMENT", 
        "\nexp:        ", exp,
        "\ncGLC:       ", cGLC,
        "\n"
    )

    epmodel_kwargs = Dict(:alpha => Inf)
    epconv_kwargs = Dict(:epsconv => 1e-4, :maxiter => 10000, :damp => 0.985)

    sim_id = string("exp: ", exp, "_global_id", SIM_GLOBAL_ID)

    get_model() = prepare_model(exp) |> ChU.uncompressed_model
    dat = ChSU.cached_simulation(;
                epochlen = EPOCHLEN, # TODO: handle better
                verbose = true,
                sim_id = sim_id,
                get_model,
                objider = iJR.BIOMASS_IDER,
                costider = iJR.COST_IDER,
                beta_info = [(iJR.BIOMASS_IDER, [β0])],
                clear_cache = false,
                use_seed = true,
                epmodel_kwargs,
                epconv_kwargs
    )

    ## CACHING RESULTS
    res_id = (:RESULTS, sim_id)
    ChU.save_cache(res_id, (exp, dat); headline = "CATCHING RESULTS\n")

    ## SAVING TO INDEX
    put!(REMCHNL, res_id)

    return nothing
end # cGLCs map

## ----------------------------------------------------------------------------
# BOUNDLING
sleep(1) # wait for RESCOLL to get all ids
const BUNDLES = Dict()
for id in RES_IDS
    exp, dat = ChU.load_cache(id; verbose = false)
    D = get!(BUNDLES, exp, Dict())
    D[:epoutb0] = dat[(:ep, 1)]
end

## ----------------------------------------------------------------------------
## SAVING
ChU.save_data(iJR.MAXENT_B0SEEDS_FILE, BUNDLES)

## ----------------------------------------------------------------------------
## CLEAR CACHE (WARNING)
if finish_clear_flag
    ChU.tagprintln_inmw("CLEARING CACHE ")
    ChU.delete_temp_caches()
end
