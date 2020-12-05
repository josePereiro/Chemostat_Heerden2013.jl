import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

## ------------------------------------------------------------------
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
## -------------------------------------------------------------------

NO_WORKERS = wcount
length(workers()) < NO_WORKERS && 
addprocs(NO_WORKERS; exeflags = "--project")
println("Working in: ", workers())

## -------------------------------------------------------------------

# Precompile in master
import Chemostat_Heerden2013

@everywhere begin
    
    import DrWatson: quickactivate 
    quickactivate(@__DIR__, "Chemostat_Heerden2013")

    using SparseArrays
    
    ## -------------------------------------------------------------------
    import Chemostat
    const ChU = Chemostat.Utils
    const ChSS = Chemostat.SteadyState
    const ChLP = Chemostat.LP
    const ChEP = Chemostat.MaxEntEP
    const ChSU = Chemostat.SimulationUtils

    ## -------------------------------------------------------------------
    import Chemostat_Heerden2013
    const Hd  = Chemostat_Heerden2013.HeerdenData;
    const BD  = Chemostat_Heerden2013.BegData;
    const iJR = Chemostat_Heerden2013.iJR904
    ChU.set_cache_dir(iJR.CACHE_DIR)

end

## ------------------------------------------------------------------
# CLEAR CACHE (WARNING)
if init_clear_flag
    ChU.tagprintln_inmw("CLEARING CACHE ")
    ChU.delete_temp_caches()
    ChU.println_inmw("\n")
end

## ------------------------------------------------------------------
# GLOBALS
@everywhere begin

    ## ------------------------------------------------------------------
    # SIMULATION GLOBAL ID
    # This must uniquely identify this simulation version
    # It is used to avoid cache collisions
    sim_global_id = "iJR904_MAXENT_EP_FBA_V2"
    
    ## -------------------------------------------------------------------
    # PARAMS
    cGLCs = Hd.val("cGLC")
    PARAMS = Dict(exp => Dict() for exp in eachindex(cGLCs))

    # This limit the beta values for each experiment,
    # for approximating the beta that enforce the experimental objective.
    # this was computed in a previous run
    for (exp, (βlb , βub)) in Dict(
                        #    βlb   ,   βub
                        1=>(56225.2, 72225.2),
                        2=>(56229.0, 72229.0),
                        3=>(49389.2, 65389.2),
                        4=>(49245.0, 65245.0),
                        5=>(61715.6, 77715.6),
                        6=>(44662.8, 60662.8),
                        7=>(44656.2, 60656.2),
                        8=>(49014.4, 65014.4),
                        9=>(49018.2, 65018.2),
                        10=>(40639.5, 56639.5),
                        11=>(40641.7, 56641.7),
                        12=>(33057.8, 49057.8),
                        13=>(33224.9, 49224.9)
                    )
        
        # Change here how many betas to model
        βcount = 30
        PARAMS[exp][:βs] = collect(range(βlb, βub, length = βcount)) 
    end

    # ep params
    # The params of EP are chosen after several tries to ensure the lower 
    # stoichiometric errors and the tuning of the objective reaction (See Cossio's paper)
    for exp in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        PARAMS[exp][:alpha] = 1e9
        PARAMS[exp][:epsconv] = 1e-6
    end

    for exp in [12,13]
        PARAMS[exp][:alpha] = 1e10
        PARAMS[exp][:epsconv] = 1e-6
    end

    # ξs
    for (exp, cGLC) in cGLCs |> enumerate
        # Change here how many xi to model, 
        # you should always include the experimental xi
        ξcount = 10
        # PARAMS[exp][:ξs] = [Hd.val("xi", exp); range(1, 100, length = ξcount)] |> sort
        PARAMS[exp][:ξs] = [Hd.val("xi", exp)] |> sort
    end

    # Intake info 
    for (exp, cGLC) in cGLCs |> enumerate
        intake_info = deepcopy(iJR.base_intake_info)

        # The only difference between experiments is the feed medium 
        # concentration.
        intake_info["EX_glc_LPAREN_e_RPAREN_"]["c"] = cGLC
            
        PARAMS[exp][:intake_info] = intake_info
    end;

    
    ## -------------------------------------------------------------------
    # PREPARE MODEL FUNCTION
    function prepare_model(ξ, intake_info)
        
        dat = ChU.load_data(iJR.BASE_MODEL_FILE; verbose = false)
        model = ChU.MetNet(;dat...)
        
        ChSS.apply_bound!(model, ξ, intake_info; 
        emptyfirst = true, ignore_miss = true)
        
        return model
    end

end # @everywhere

## -------------------------------------------------------------------
# SIM IDS
# Collect all the computed simulation ids for collecting results
const chnl = RemoteChannel() do
    Channel{Any}(10)
end
const res_ids = []
const collector = @async while true
    id = take!(chnl)
    push!(res_ids, id)
end

## -------------------------------------------------------------------
# RUN SIMULATION
empty!(res_ids)
pmap(enumerate(cGLCs)) do (exp, cGLC) # This is parallizable

    EPOCHLEN = 100

    exp_params = PARAMS[exp]
    ξs = exp_params[:ξs]
    βs = exp_params[:βs]
    intake_info = exp_params[:intake_info]

    epmodel_kwargs = Dict(k => exp_params[k] for k in [:alpha])
    epconv_kwargs = Dict(k => exp_params[k] for k in [:epsconv])

    ChU.tagprintln_inmw("PROCESSING EXPERIMENT", 
        "\nexp:        ", exp,
        "\ncGLC:       ", cGLC,
        "\nxi count:   ", length(ξs),
        "\nbeta count: ", length(βs),
        "\nalpha:      ", exp_params[:alpha],
        "\nepsconv:    ", exp_params[:epsconv],
        "\nEPOCHLEN:   ", EPOCHLEN,
        "\n"
    )


    for (ξi, ξ) in ξs |> enumerate

        sim_id = string("exp: ", exp, " xi: [", ξi, "/", length(ξs), "] ",
                    "global_id: ", sim_global_id, "_", hash(PARAMS))

        ChU.tagprintln_inmw("PROCESSING XI", 
            "\nxi: ", ξ, " [", ξi, ",", length(ξs), "]",
            "\n"
        )

        dat = ChSU.cached_simulation(;
                    epochlen = EPOCHLEN, # TODO: handle better
                    verbose = true,
                    sim_id = sim_id,
                    get_model = () -> prepare_model(ξ, intake_info),
                    objider = iJR.BIOMASS_IDER,
                    costider = iJR.COST_IDER,
                    beta_info = [(iJR.BIOMASS_IDER, βs)],
                    clear_cache = false,
                    use_seed = true,
                    epmodel_kwargs = epmodel_kwargs,
                    epconv_kwargs = epconv_kwargs
        )

        ## CACHING RESULTS
        res_id = (:RESULTS, sim_id)
        model = prepare_model(ξ, intake_info)
        ChU.save_cache(res_id, (exp, ξ, βs, model, dat); 
            headline = "CATCHING RESULTS\n")

        ## SAVING TO INDEX
        put!(chnl, res_id)

    end # ξs

    return nothing
end # cGLCs map

## -------------------------------------------------------------------
# BOUNDLING
sleep(1) # wait for collector to get all ids
bundles = Dict()
for id in res_ids
    exp, ξ, βs, model, dat = ChU.load_cache(id; verbose = false)

    # boundling
    bundle = get!(bundles, exp, ChU.ChstatBundle())

    bundle[ξ, :net] = model
    bundle[ξ, :fba] = dat[:fba]

    for (βi, β) in βs |> enumerate
        bundle[ξ, β, :ep] = dat[(:ep, βi)]
    end
end

## SAVING
ChU.save_data(iJR.MAXENT_FBA_EB_BOUNDLES_FILE, bundles)

## CLEAR CACHE (WARNING)
if finish_clear_flag
    ChU.tagprintln_inmw("CLEARING CACHE ")
    ChU.delete_temp_caches()
end
