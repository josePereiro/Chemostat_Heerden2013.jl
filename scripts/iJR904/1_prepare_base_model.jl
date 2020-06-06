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

import CSV
using DataFrames
using Serialization

# run add "https://github.com/josePereiro/Chemostat" in the 
# julia Pkg REPL for installing the package
import Chemostat
Ch = Chemostat

# +
# Run add https://github.com/josePereiro/Chemostat_Heerden2013.jl in the Julia Pkg REPL to install the
# package, then you must activate the package enviroment (see README)
import Chemostat_Heerden2013
Hd  = Chemostat_Heerden2013.HeerdenData;
Bd  = Chemostat_Heerden2013.BegData
iJR = Chemostat_Heerden2013.iJR904

# This just check that the script is run in the
# package enviroment
Chemostat_Heerden2013.check_env()
# -

# ---
# ## BASE MODEL
# ---

# ### Load Mat file

# +
println("Original .mat model")
model = Ch.Utils.read_mat(iJR.MODEL_RAW_MAT_FILE);
obj_ider = "BiomassEcoli"

# the maximal experimental growth rate in Heerden2013 is ~0.2 1/h
# The raw model present a growth rate bigger than that, so it is ok
# to use it directly as base model
fbaout = Ch.FBA.fba(model, obj_ider);
Ch.Utils.summary(model, fbaout)
# -

# ### Set bounds

# +
## The abs maximum bounds will be set to 100
abs_max_bound = 100

foreach(function (ider)
        Ch.Utils.isfixxed(model, ider) && return # fixxed reaction are untouched
        old_ub = Ch.Utils.ub(model, ider)
        new_ub = old_ub == 0.0 ? 0.0 : abs_max_bound
        Ch.Utils.ub!(model, ider, new_ub)
    end, model.rxns);

foreach(function (ider)
        Ch.Utils.isfixxed(model, ider) && return # fixxed reaction are untouched
        old_lb = Ch.Utils.lb(model, ider)
        new_lb = old_lb == 0.0 ? 0.0 : -abs_max_bound
        Ch.Utils.lb!(model, ider, new_lb)
        end, model.rxns);

Ch.Utils.summary(model)
# -

# ### Enzymatic cost info

# +
# Close, for now, all exchanges for avoiding it to be in revs
# The reversible reactions will be splited for modeling cost
# Exchanges have not associated cost, so, we do not split them
foreach(idx -> Ch.Utils.ub!(model, idx, 0.0), 
    Ch.Utils.exchanges(model)) # Closing all outtakes
foreach(idx -> Ch.Utils.lb!(model, idx, 0.0), 
    Ch.Utils.exchanges(model)) # Closing all intakes


# The cost will be introduced as a reaction, we follow the same cost models as 
# Beg et al. (2007): https://doi.org/10.1073/pnas.0609845104.
# A new balance equations is then added:
#        Σ(rᵢ*costᵢ) + tot_cost = 0
#    Because the cost coefficients (costᵢ) < 0 (it resamble a reactant), the system must allocate 
#    the fluxes (rᵢ) so that Σ(rᵢ*costᵢ) = tot_cost, and tot_cost
#    are usually bounded [0.0, 1.0]
cost_info = Dict()
for rxn in model.rxns
    # The exchanges, the atpm and the biomass are synthetic reactions, so, 
    # they have should not have an associated enzimatic cost 
    any(startswith.(rxn, ["EX_", "DM_"])) && continue
    rxn == "BiomassEcoli" && continue
    rxn == "ATPM" && continue
        
    # We will split the rev reactions, so we save the cost for both versions (fwd, bkwd)
    if Ch.Utils.isrev(model, rxn)
        cost_info[rxn * "_fwd"] = -iJR.beg_enz_cost(rxn)
        cost_info[rxn * "_bkwd"] = -iJR.beg_enz_cost(rxn)
        continue
    end
    
    
    # Only the internal, non reversible reactions have an associated cost
    cost_info[rxn] = -iJR.beg_enz_cost(rxn)
end

# Spliting revs
model = Ch.Utils.split_revs(model);

# Adding cost raction
model = Ch.Utils.add_costs(model, cost_info, verbose = false);
# -

# ### Set base exchanges

# To control the intakes just the metabolites defined in the 
# base_intake_info (The minimum medium) will be opened.
# The base model will be constraint as in a cultivation with xi = 1.0
# see Cossios paper (see README)
ξ = 1.0 
println("Minimum medium: ", iJR.base_intake_info)
foreach(idx -> Ch.Utils.lb!(model, idx, 0.0), 
    Ch.Utils.exchanges(model)) # Closing all intakes
foreach(idx -> Ch.Utils.ub!(model, idx, 100.0), 
    Ch.Utils.exchanges(model)) # Opening all outakes all
# see Cossios paper (see README) for details in the Chemostat bound constraint
Ch.SteadyState.apply_bound!(model, ξ, iJR.base_intake_info);

# tot_cost is the exchange that controls the bounds of the 
# enzimatic cost contraint, we bound it to [0, 1.0]
Ch.Utils.lb!(model, "tot_cost", 0.0);
Ch.Utils.ub!(model, "tot_cost", 1.0);

# ### TODO: Delete, this is just for checking

# Checking with old gem
old_rxns_path = "/Users/Pereiro/University/Physic/Research/Jose/Metabolism/MaxEntEP/projects/EColi_Heerden2013/data/processed/gems/iJR904___rxns.tsv"
old_rxns = DataFrame(CSV.read(old_rxns_path, delim = "\t"));
old_mets_path = "/Users/Pereiro/University/Physic/Research/Jose/Metabolism/MaxEntEP/projects/EColi_Heerden2013/data/processed/gems/iJR904___mets.tsv"
old_mets = DataFrame(CSV.read(old_mets_path, delim = "\t"));


# +
# Cheking internal reactions
cost_met = "cost"
for old_rxn in eachrow(old_rxns)
    
    # We must consider the split rev rename
    new_idxs = [findfirst(isequal(old_rxn.id), model.rxns), 
                findfirst(isequal(old_rxn.id * "_fwd"), model.rxns),
                findfirst(isequal(old_rxn.id * "_bkwd"), model.rxns)]
    
    all(isnothing.(new_idxs)) && (println("$(old_rxn.id) missing"); continue)

    for (i, new_idx) in enumerate(new_idxs)
        
        isnothing(new_idx) && continue
        
        new_id = model.rxns[new_idx]
        # Bounds
        new_lb, new_ub = Ch.Utils.bounds(model, new_idx)
        if i == 1 # Non rev
            old_lb, old_ub = old_rxn.nub, old_rxn.pub
        elseif i == 2 # fwd
            old_lb, old_ub = 0.0, old_rxn.pub
            else # bkwd
            old_lb, old_ub = 0.0, -old_rxn.nub
        end

        !(new_lb ≈ old_lb) && println("$new_id: not approx lb new: ", 
            new_lb, " old: ", old_lb)
        !(new_ub ≈ old_ub) && println("$new_id: not approx ub new: ", 
            new_ub, " old: ", old_ub)
        !(new_ub ≈ old_ub) && println()

        # Cost
        cost_idx = findfirst(isequal(cost_met), model.mets)
        isnothing(new_idx) && continue

        old_rxn.ap != old_rxn.an &&
            println(old_rxn.ap, " different ap, an", old_rxn.ap ,old_rxn.an)

        if cost_idx in Ch.Utils.rxn_reacts(model, new_idx)
            new_a = Ch.Utils.S(model, cost_idx, new_idx)
            old_a = old_rxn.ap
            if !(abs(new_a) ≈ abs(old_a))
                println("$new_id: not approx a new: ", 
            new_a, " old: ", old_a)
            end
        else
            println("new $new_id have not defined cost, old ap: ", old_rxn.ap, " old an: ", old_rxn.an)
        end
    end
    
end
# -

### Biomass
obj_ider = "BiomassEcoli"
for new_idx in Ch.Utils.rxn_mets(model, obj_ider)
    new_met = model.mets[new_idx]
    old_idx = findfirst(isequal(new_met), old_mets.id)
    new_s = Ch.Utils.S(model, new_idx, obj_ider)
    old_s = old_mets[old_idx, :y]
    !(old_s ≈ -new_s) && println("$new_met: not approx y new: ", 
        new_s, " old: ", old_s)
end

### ATPM
ider = "ATPM"
for new_idx in Ch.Utils.rxn_mets(model, ider)
    new_met = model.mets[new_idx]
    old_idx = findfirst(isequal(new_met), old_mets.id)
    new_s = Ch.Utils.S(model, new_idx, obj_ider)
    old_s = old_mets[old_idx, :y]
    !(old_s ≈ -new_s) && println("$new_met: not approx y new: ", 
        new_s, " old: ", old_s)
end

# Cheking external reactions
for old_met in eachrow(old_mets)
    # Medium
    if old_met.L != 0.0 || old_met.V != 0.0 
        new_exch = get(iJR.exch_met_map, old_met.id, nothing)
        isnothing(new_exch) && (println("$(old_met.id) not found"); continue)
        lb, ub = Ch.Utils.bounds(model, new_exch)
        L, V = old_met.L, old_met.V
        if !(abs(lb) ≈ abs(V))
            println("new $new_exch, lb: $lb")
            println("old $(old_met.id), V: $V")
            println()
        end
        if !(abs(ub) ≈ abs(L))
            println("new $new_exch, ub: $ub")
            println("old $(old_met.id), L: $L")
            println()
        end
    end
end

# ---
# ## Exch_met_map
# ---

# +
# A quick way to que exchages from mets and te other way around
exch_met_map = Dict()
for exch_ in model.rxns[findall((id) -> any(startswith.(id, ["EX_", "DM_"])), model.rxns)]
    mets_ = model.mets[Ch.Utils.rxn_mets(model, exch_)]
    length(mets_) != 1 && continue
    exch_met_map[exch_] = mets_[1]
    exch_met_map[mets_[1]] = exch_
end;

# saving
df = DataFrame(x1 = collect(keys(exch_met_map)), x2 = collect(values(exch_met_map)))
CSV.write(iJR.EXCH_MET_MAP_FILE, df)
println(relpath(iJR.EXCH_MET_MAP_FILE), " created!!!")
# -

# ### Saving model

# +
model = Ch.Utils.fva_preprocess(model, 
#     eps = 1-9, # This avoid blocking totally any reaction
    verbose = true);

fbaout = Ch.FBA.fba(model, obj_ider);
Ch.Utils.summary(model, fbaout)
Ch.Utils.summary(model)

## Saving model
serialize(iJR.BASE_MODEL_FILE, model)
println(relpath(iJR.BASE_MODEL_FILE), " created!!!")
