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

NO_CORES = 3#length(Sys.cpu_info())
# length(workers()) < NO_CORES && addprocs(NO_CORES)
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

notebook_name = "maxent_ep_v1";

# ### Params

# +
@everywhere begin
params = Dict()

# Intake info
params["obj_ider"] = "BiomassEcoli";
params["cost_ider"] = "tot_cost";

# Intake info
params["intake_infos"] = []
for (i, cGLC) in enumerate(Hd.val("cGLC"))
    intake_info = deepcopy(iJR.base_intake_info)
    intake_info["EX_glc_LPAREN_e_RPAREN_"]["c"] = cGLC
        
    push!(params["intake_infos"], intake_info)
end;

# ep params
params["ep_epsconv"] = 1e-6
params["ep_α"] = 1e9

# This limit the beta values for each experiment,
# for approximating the beta that enforce the experimental objective.
# this was computed in a previous run
params["ep_β_intervals"] = Dict(
    1=>(11724.5, 16724.5),
    2=>(11724.5, 16724.5),
    3=>(28153.1, 33153.1),
    4=>(19704.1, 24704.1),
    5=>(14071.4, 19071.4),
    7=>(7500.0, 12500.0),
    6=>(7500.0, 12500.0),
    8=>(7969.39, 12969.4),
    9=>(7969.39, 12969.4),
    10=>(8908.16, 13908.2),
    11=>(8908.16, 13908.2),
    12=>(27214.3, 32214.3),
    13=>(30500.0, 35500.0)
)

    
end
# -

# ### Print functions

@everywhere function print_hello(wid, exp, ξs, βs, t)
    Core.println("Worker $wid starting exp $exp at $(t) ----------------------------")
    Core.println("xis: $ξs")
    Core.println("betas: $βs")
    Core.println()
    flush(stdout);
end

@everywhere function print_progress(wid, exp, ξi, ξs, ξ,  βi, βs, β, exp_av, fba_av, ep_av, elapsed)
    Core.println("w: $wid exp: $exp progress ----------------------------")
    Core.println("\texp cGLC: $(Hd.val(:cGLC, exp))")
    Core.println("\texp xi: $(Hd.val(:xi, exp))")
    Core.println("\txi  [$ξi, $(length(ξs))]: $ξ")
    Core.println("\tbeta[$βi, $(length(βs))]: $β")
    Core.println("\texp_av: $exp_av")
    Core.println("\tfba_av: $fba_av")
    Core.println("\tep_av : $ep_av")
    Core.println("\ttime: $elapsed second(s)")
    Core.println()
    flush(stdout);
end

@everywhere function print_good_bye(wid, exp, t)
    Core.println("Worker $wid finished exp $exp at $(t)----------------------------")
    Core.println()
    flush(stdout);
end

# ### work_function

@everywhere function process_exp(exp)

    # ep params
#     ξs = [Hd.val("xi", exp); round.(Ch.Utils.logspace(-1, 1, 100), digits = 3)] |> sort
    ξs = Ch.Utils.logspace(-1, 2, 20)
#     βlb, βub = params["ep_β_intervals"][exp]
#     βs = collect(range(βlb, βub, length = 50))
    βs = []
    
    cGLC = round(Hd.val("cGLC", exp), digits = 3)
    
    # Info
    # Print hello in worker 1
#     remotecall_wait(print_hello, 1, myid(), exp, ξs, βs, now())
    
    boundle = Ch.Utils.ChstatBoundle()
    
    for (ξi, ξ) in ξs |> sort |> reverse |> enumerate 
        # Prepare model
        model = deserialize(iJR.BASE_MODEL_FILE);
        
        # chsts bound
        Ch.SteadyState.apply_bound!(model, ξ, params["intake_infos"][exp])
        obj_idx = Ch.Utils.rxnindex(model, params["obj_ider"])
        
        # fba
        fbaout = Ch.FBA.fba(model, params["obj_ider"], params["cost_ider"])
#         fbaout = Ch.FBA.fba(model, params["obj_ider"])
        
        
        # boundling
        Ch.Utils.add_data!(boundle, ξ, :fba, fbaout)
        Ch.Utils.add_data!(boundle, ξ, :net, model)
        
        # Info
        exp_av = Hd.val("D", exp)
        fba_av = Ch.Utils.av(model, fbaout, obj_idx)
        
        t0 = time() # to track duration
        seed_epout = nothing
        βv = zeros(size(model, 2))
        for (βi, β) in βs |> sort |> reverse |> enumerate
            try
                # epout
                βv[obj_idx] = β
                epout = Ch.MaxEntEP.maxent_ep(model; 
                                            α = params["ep_α"], 
                                            βv = βv,
                                            epsconv = params["ep_epsconv"], 
                                            solution = seed_epout,
                                            verbose = false)
                seed_epout = epout
                
                # Info
                ep_av = Ch.Utils.av(model, epout, obj_idx)
                
                # Print progress in worker 1
                βi % 10 == 0 && remotecall_wait(print_progress, 1, myid(), exp, ξi, ξs, ξ,  βi, βs, 
                    β, exp_av, fba_av, ep_av, time() - t0)
        
                Ch.Utils.add_data!(boundle, ξ, β, :ep, epout)
                
            catch err
                println("Error at exp: $exp, xi: $ξ, beta; $β")
                rethrow(err)
            end
        end
        
    end
    
    # Info
    # Print hello in worker 1
    remotecall_wait(print_good_bye, 1, myid(), exp, now())
    
    return (exp, boundle)
end

# ---
# ### Parallel loop
# ---
# This will output a lot of test and can take a while!!!

# remote_results = pmap(process_exp, eachindex(Hd.val("cGLC")));
remote_results = map(process_exp, eachindex(Hd.val("cGLC")));

boundles = Dict();
for (exp, boundle) in remote_results
    boundles[exp] = boundle
end

# ---
# ### Saving
# ---

# +
# cache_file = joinpath(iJR.MODEL_PROCESSED_DATA_DIR, "$(notebook_name)__cache2.jls")
# serialize(cache_file, (params, boundles))
# println(realpath(cache_file), " created!!!")
# -
# ---
# ## Checking
# ---

using Plots

colors = Plots.distinguishable_colors(length(boundles))

# +
# TODO package this
# Transforme a exchange value to a medium concentration (mM), it is by using the 
# steady state assumption in the chemostat, see Cossio's paper (see README)
function conc(Hd_met, exp, ξ::Real, data_idxs...)
    # s = c - u*ξ if u > 0 means intake
    model_met = iJR.Hd_mets_map[Hd_met]
    model_exch = iJR.exch_met_map[model_met]
    boundle = boundles[exp]
    
    u = Ch.Utils.av(boundle, ξ, data_idxs..., model_exch) # exchange val
    uerr = Ch.Utils.va(boundle, ξ, data_idxs..., model_exch) |> sqrt # exchange std
    intake_info = params["intake_infos"][exp]
    c = Hd_met == "GLC" ? Hd.val("cGLC", exp) : 0.0
    return (max(c + u*ξ, 0.0), uerr * ξ)
end

function conc(Hd_met, exp, ξs::Vector, data_idxs...)
    model_concs, model_conc_stds = [], []
    foreach(function (ξ) 
                model_conc, model_conc_std = conc(Hd_met, exp, ξ, data_idxs...)
                push!(model_concs, model_conc)
                push!(model_conc_stds, model_conc_std)
            end, ξs)
    return model_concs, model_conc_stds
end

# +
# exp = 1
# @show ξ = boundles[exp].ξs[1]
# model = Ch.Utils.get_data(boundles[exp], ξ, :net);
# Ch.Utils.summary(model)
# -

# Biomass
ider = params["obj_ider"]
p = Plots.plot(title = ider, xlabel = "xi", ylabel = "flx", 
#     xaxis = [20, 25], 
#     yaxis = [0.4, 0.5]
)
for (exp, boundle) in boundles
    avs = Ch.Utils.av(boundle, boundle.ξs, :fba, ider)
    Plots.plot!(p, boundle.ξs, avs, color = colors[exp], label = "")
    Plots.scatter!(p, [Hd.val("xi", exp)], [Hd.val("D", exp)], color = colors[exp], label = "")
end
p

# Reactions
ider = Hd.msd_mets[1]
ider = iJR.exch_met_map[iJR.Hd_mets_map[ider]] # * "_bkwd"
p = Plots.plot(title = ider, xlabel = "xi", ylabel = "flx", 
#     yaxis = [0.0, 25.0],
#         xaxis = [25, 30]
)
for (exp, boundle) in boundles
    avs = -Ch.Utils.av(boundle, boundle.ξs, :fba, ider)
    Plots.plot!(p, boundle.ξs, avs, color = colors[exp], label = "")
    models = Ch.Utils.get_data(boundle, boundle.ξs, :net)
    lbs = -map(model -> Ch.Utils.lb(model, ider), models)
#     Plots.plot!(p, boundle.ξs, lbs, ls = :dash,  color = colors[exp], label = "")
end
p

model = deserialize(iJR.BASE_MODEL_FILE);
Ch.Utils.summary(model, "ATPM")

# Cost
cost_met = "cost"
cost_rxn = "tot_cost"
p = Plots.plot(title = cost_rxn, xlabel = "xi", ylabel = "cost",
    yaxis = [0.0, 1.0]
)
for (exp, boundle) in boundles
    avs = Ch.Utils.av(boundle, boundle.ξs, :fba, cost_rxn)
    Plots.plot!(p, boundle.ξs, avs, color = colors[exp], label = "")
end
p



# Metabolites Concentration
Hd_ider = Hd.msd_mets[1]
p = Plots.plot(title = Hd_ider, xlabel = "xi", ylabel = "conc",
#         xaxis = [0.0, 5.0],
    );
for (exp, boundle) in boundles
    concs_, conc_errs_ = conc(Hd_ider, exp, boundle.ξs, :fba)
    Plots.plot!(p, boundle.ξs, concs_, color = colors[exp], label = "")
    Plots.scatter!(p, [Hd.val("xi", exp)], [Hd.val("s$Hd_ider", exp)], color = colors[exp], label = "")
end
p


