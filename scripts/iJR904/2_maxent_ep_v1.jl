
using Distributed

NO_CORES = length(Sys.cpu_info()) - 1
length(workers()) < NO_CORES && addprocs(NO_CORES - 1; exeflags = "--project")
println("Working in: ", workers())

@everywhere begin
    
    using DrWatson 
    quickactivate(@__DIR__, "Chemostat_Heerden2013")

    using SparseArrays
    
    ## -------------------------------------------------------------------
    import Chemostat
    import Chemostat.Utils: MetNet, EPModel, ChstatBoundle,
                            rxnindex, metindex, compress_dict, 
                            uncompress_dict, clampfileds!, well_scaled_model,
                            ChstatBoundle, norm_abs_stoi_err, av, va, nzabs_range,
                            struct_to_dict
    import Chemostat.SimulationUtils: epoch_converge_ep!, cached_simulation, set_cache_dir, 
                            tagprintln_inmw, println_inmw, tagprintln_ifmw, println_ifmw,
                            save_cache, load_cache, delete_temp_caches
    import Chemostat.SteadyState: apply_bound!

    ## -------------------------------------------------------------------
    import Chemostat_Heerden2013: HeerdenData, BegData, iJR904, save_data, load_data
    import Chemostat_Heerden2013.iJR904: OBJ_IDER, ATPM_IDER, COST_IDER
    const Hd  = HeerdenData;
    const iJR = iJR904
    set_cache_dir(iJR.MODEL_PROCESSED_DATA_DIR)

end


@everywhere begin

    ## ------------------------------------------------------------------
    # SIMULATION GLOBAL ID
    # This must uniquely identify this simulation version
    # It is used to avoid cache collisions
    sim_global_id = "iJR904_maxent_fba_ep_v1"
    
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
        PARAMS[exp][:βs] = collect(range(βlb, βub, length = 30)) 
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
        # Change here how many xi to model, you should always include the experimental xi
        PARAMS[exp][:ξs] = [Hd.val("xi", exp); range(1, 100, length = 10)] |> sort
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
        
        dat = load_data(iJR.BASE_MODEL_FILE; verbose = false)
        model = MetNet(;dat...)
        
        apply_bound!(model, ξ, intake_info; 
        emptyfirst = true, ignore_miss = true)
        
        return model
    end

end # @everywhere

## -------------------------------------------------------------------
# RUN SIMULATION

windexes = pmap(enumerate(cGLCs)) do (exp, cGLC) # This is parallizable

    EPOCHLEN = 40
    index_id = (:INDEX, myid())

    exp_params = PARAMS[exp]
    ξs = exp_params[:ξs]
    βs = exp_params[:βs]
    intake_info = exp_params[:intake_info]

    epmodel_kwargs = Dict(k => exp_params[k] for k in [:alpha])
    epconv_kwargs = Dict(k => exp_params[k] for k in [:epsconv])

    tagprintln_inmw("PROCESSING EXPERIMENT", 
        "\nexp:        ", exp,
        "\ncGLC:       ", cGLC,
        "\nxi count:   ", length(ξs),
        "\nbeta count: ", length(βs),
        "\nalpha:      ", exp_params[:alpha],
        "\nepsconv:    ", exp_params[:epsconv],
        "\nEPOCHLEN:   ", EPOCHLEN,
        "\n"
    )

    exp_id = (exp, cGLC, sim_global_id, hash(PARAMS))

    for (ξi, ξ) in ξs |> enumerate

        sim_id = (ξi, ξ, exp_id)

        tagprintln_inmw("PROCESSING XI", 
            "\nxi: ", ξ, " [", ξi, ",", length(ξs), "]",
            "\n"
        )

        dat = cached_simulation(;
                    epochlen = EPOCHLEN, # TODO: handle better
                    verbose = true,
                    sim_id = sim_id,
                    get_model = function()
                        return prepare_model(ξ, intake_info);
                    end,
                    objider = OBJ_IDER,
                    costider = COST_IDER,
                    beta_info = [(OBJ_IDER, βs)], # TODO: handle betas better 
                    clear_cache = false,
                    use_seed = true,
                    epmodel_kwargs = epmodel_kwargs,
                    epconv_kwargs = epconv_kwargs
        )

        ## CACHING RESULTS
        model = prepare_model(ξ, intake_info)
        save_cache(sim_id, (exp, ξ, βs, model, dat); headline = "CATCHING RESULTS\n")

        ## SAVING TO INDEX
        index = load_cache(index_id; verbose = false)
        index = isnothing(index) ? [] : index
        push!(index, sim_id)
        save_cache(index_id, index; verbose = false)

    end # ξs

    return index_id
end # cGLCs map


## BOUNDLING
boundles = Dict()
for index_id in windexes

    index = load_cache(index_id; verbose = false)
    isnothing(index) && continue
    
    for id in index
        exp, ξ, βs, model, dat = load_cache(id; verbose = false)
       
        # boundling
        boundle = get!(boundles, exp, ChstatBoundle())

        boundle[ξ, :net] = model
        boundle[ξ, :fba] = dat[:fba]

        for (βi, β) in βs |> enumerate
            boundle[ξ, β, :ep] = dat[(:ep, βi)]
        end

    end
end

## SAVING
save_data(iJR.MODEL_PROCESSED_DATA_DIR, boundles)

## DELETING CACHE
delete_temp_caches()

