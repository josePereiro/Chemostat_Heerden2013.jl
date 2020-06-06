# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: jl,ipynb
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Julia 1.1.0
#     language: julia
#     name: julia-1.1
# ---

# ### Precompaling in master worker first

# run add "https://github.com/josePereiro/Chemostat" in the 
# julia Pkg REPL for installing the package
import Chemostat
Ch = Chemostat

# +
# Run add https://github.com/josePereiro/Chemostat_Heerden2013.jl in the Julia Pkg REPL to install the
# package, then you must activate the package enviroment (see README)
import Chemostat_Heerden2013
Hd = Chemostat_Heerden2013.HeerdenData;
Bd = Chemostat_Heerden2013.BegData
iJR = Chemostat_Heerden2013.iJR904

# This just check that the script is run in the
# package enviroment
Chemostat_Heerden2013.check_env()
# -

# ### Loading in everywhere

# +
using Distributed

NO_CORES = length(Sys.cpu_info())
length(workers()) < NO_CORES && addprocs(NO_CORES - 1)
println("Working in: ", workers())

# +
@everywhere begin

using Distributed
using Serialization
using Dates

using Chemostat
Ch = Chemostat
    
# Run add https://github.com/josePereiro/Chemostat_Heerden2013.jl in the Julia Pkg REPL to install the
# package, then you must activate the package enviroment (see README)
import Chemostat_Heerden2013
Hd = Chemostat_Heerden2013.HeerdenData;
Bd = Chemostat_Heerden2013.BegData
iJR = Chemostat_Heerden2013.iJR904

# This just check that the script is run in the
# package enviroment
Chemostat_Heerden2013.check_env()
    
end
# -

#
# ### Meta

@everywhere notebook_name = "maxent_ep_v1";

# ### Params

# +
@everywhere begin

params = Dict()

# Intake info
params["obj_ider"] = "BiomassEcoli";
params["cost_ider"] = "tot_cost";

# Intake info
# TODO redesign this to be included in the model
params["intake_infos"] = Dict()
for (exp, cGLC) in enumerate(Hd.val("cGLC"))
    intake_info = deepcopy(iJR.base_intake_info)

    # The only difference between experiments is the feed medium 
    # concentration.
    intake_info["EX_glc_LPAREN_e_RPAREN_"]["c"] = cGLC
        
    params["intake_infos"][exp] = intake_info
end;

# ep params
# The params of EP are chosen after several tries to ensure the lower 
# stoichiometric errors and the tuning of the objective reaction (See Cossio's paper)
params["ep_alpha"] = Dict()
params["ep_epsconv"] = Dict()
for exp in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    params["ep_alpha"][exp] = 1e9
    params["ep_epsconv"][exp] = 1e-6
end
for exp in [12,13]
    params["ep_alpha"][exp] = 1e10
    params["ep_epsconv"][exp] = 1e-6
end

# This limit the beta values for each experiment,
# for approximating the beta that enforce the experimental objective.
# this was computed in a previous run
# TODO automatize this
params["ep_beta_intervals"] = Dict(
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
 
end
# -

# ### Print functions

@everywhere function print_hello(wid, exp, ξs, βs)
    Core.println("Worker $wid starting exp $exp at $(Time(now())) ----------------------------")
    Core.println("\txis:   $ξs")
    Core.println("\tbetas: $βs")
    Core.println()
    flush(stdout);
end

@everywhere function print_progress(wid, exp, ξi, ξs, ξ,  βi, βs, β, 
            exp_av, fba_av, ep_av, ep_alpha, ep_epsconv,
            elapsed)
    Core.println("Worker: $wid exp: $exp progress at $(Time(now())) ----------------------------")
    Core.println("\txi: [$ξi/ $(length(ξs))] beta: [$βi/ $(length(βs))]")
    Core.println("\t  ----------------- --------------")
    Core.println("\tmodel xi:           $ξ")
    Core.println("\tmodel beta:         $β")
    Core.println("\tmodel ep_alpha:     $ep_alpha")
    Core.println("\tmodel ep_epsconv:   $ep_epsconv")
    Core.println("\tmodel fba obj:      $fba_av")
    Core.println("\tmodel ep obj:       $ep_av")
    Core.println("\t  ----------------- --------------")
    Core.println("\texp   xi:           $(Hd.val(:xi, exp))")
    Core.println("\texp   cGLC:         $(Hd.val(:cGLC, exp))")
    Core.println("\texp   obj:          $exp_av")
    Core.println("\t  ----------------- --------------")
    Core.println("\txi elapsed time(s): $elapsed")
    Core.println()
    flush(stdout);
end

@everywhere function print_good_bye(wid, exp, tcache_file)
    Core.println("Worker $wid finished exp $exp at $(Time(now())) ----------------------------")
    Core.println("\tTemp cache file ", relpath(tcache_file))
    Core.println()
    flush(stdout);
end

@everywhere function print_return_cached(wid, exp, tcache_file)
    Core.println("Worker $wid returns cached exp $exp at $(Time(now())) ----------------------------")
    Core.println("\tTemp cache file ", relpath(tcache_file))
    Core.println()
    flush(stdout);
end

# ### temp cache file

@everywhere temp_cache_file(exp) = 
    joinpath(iJR.MODEL_PROCESSED_DATA_DIR, "$(notebook_name)_temp_cache_exp$exp.jls")

# ### work_function

@everywhere function process_exp(exp; upfrec = 10)

    # I will cache temporally the results of this function 
    # I do not check any cache consistency, so delete the temporal caches if you
    # dont trust them
    tcache_file = temp_cache_file(exp)

    # If cached return 
    if isfile(tcache_file)
        # Info
        # Print in worker 1
        remotecall_wait(print_return_cached, 1, myid(), exp, tcache_file)

        return deserialize(tcache_file)
    end

    # prepare params
    # Change here how many xi to model, you should always include the experimental xi
    ξs = [Hd.val("xi", exp); range(1, 100, length = 10)] |> sort 

    βlb, βub = params["ep_beta_intervals"][exp]
    # Change here how many betas to model
    βs = collect(range(βlb, βub, length = 30)) 

    ep_alpha = params["ep_alpha"][exp]
    ep_epsconv = params["ep_epsconv"][exp]

    obj_ider = params["obj_ider"]
    cost_ider = params["cost_ider"]
    intake_info = params["intake_infos"][exp]
    
    # Info
    # Print hello in worker 1
    remotecall_wait(print_hello, 1, myid(), exp, ξs, βs)
    
    boundle = Ch.Utils.ChstatBoundle()
    
    for (ξi, ξ) in ξs |> sort |> reverse |> enumerate 
        # Prepare model
        model = deserialize(iJR.BASE_MODEL_FILE);

        # chsts bound
        Ch.SteadyState.apply_bound!(model, ξ, intake_info)
        obj_idx = Ch.Utils.rxnindex(model, obj_ider)
        
        # fba
        fbaout = Ch.FBA.fba(model, obj_ider, cost_ider)
        
        # boundling
        Ch.Utils.add_data!(boundle, ξ, :fba, fbaout)
        Ch.Utils.add_data!(boundle, ξ, :net, model)
        
        # Info
        exp_av = Hd.val("D", exp)
        fba_av = Ch.Utils.av(model, fbaout, obj_idx)
        
        t0 = time() # to track xi processing duration
        seed_epout = nothing
        βv = zeros(size(model, 2))
        for (βi, β) in βs |> sort |> reverse |> enumerate
            try
                # epout
                βv[obj_idx] = β
                epout = Ch.MaxEntEP.maxent_ep(model; 
                                            α = ep_alpha, 
                                            βv = βv,
                                            epsconv = ep_epsconv, 
                                            solution = seed_epout,
                                            verbose = false)
                # The previous result is used as starting point in the next computation, 
                # this reduce the convergence time a lot 
                seed_epout = epout 
                
                

                # Boundling
                Ch.Utils.add_data!(boundle, ξ, β, :ep, epout)
                
                # Info
                # Print progress in worker 1
                ep_av = Ch.Utils.av(model, epout, obj_idx)
                show_progress = βi == 1 || βi == length(βs) || βi % upfrec == 0
                show_progress && remotecall_wait(print_progress, 1, myid(), 
                    exp, ξi, ξs, ξ,  βi, βs, β, 
                    exp_av, fba_av, ep_av, ep_alpha, ep_epsconv,
                    time() - t0)
                
            catch err
                println("Worker $(myid()): Error at exp: $exp, xi: $ξ, beta; $β")
                rethrow(err)
            end
        end
        
    end
    
    # Info
    # Print hello in worker 1
    remotecall_wait(print_good_bye, 1, myid(), exp, tcache_file)
    
    # Catching
    serialize(tcache_file, (exp, boundle))

    return (exp, boundle)
end

# ---
# ### Parallel loop
# ---
# This will output a lot of text and can take a while (~ 4h, 3 workers, Dual-Core Intel Core i5)!!!

remote_results = pmap(process_exp, eachindex(Hd.val("cGLC")));


boundles = Dict();
for (exp, boundle) in remote_results
    boundles[exp] = boundle
end

# ---
# ### Saving results cache
# ---

println("Saving results")
cache_file = joinpath(iJR.MODEL_PROCESSED_DATA_DIR, "$(notebook_name)__cache8.jls")
serialize(cache_file, (params, boundles))
println(realpath(cache_file), " created!!!")
println()

# ---
# ### Deleting temp caches
# ---

println("Deleting temporal cache files")
for exp in eachindex(Hd.val("cGLC"))
    tcache_file = temp_cache_file(exp)
    rm(tcache_file, force = true)
    println(relpath(tcache_file), " deleted!!!")
end
