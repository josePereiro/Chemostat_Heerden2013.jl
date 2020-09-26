using DrWatson 
quickactivate(@__DIR__, "Chemostat_Heerden2013")

import CSV
import MAT
using DataFrames
using Serialization

## -------------------------------------------------------------------
# run add "https://github.com/josePereiro/Chemostat" in the 
# julia Pkg REPL for installing the package
import Chemostat
import Chemostat.Utils: MetNet, to_symbol_dict, isrev, split_revs, rxn_mets,
                        rxnindex, metindex, compress_dict, exchanges, expanded_model,
                        uncompress_dict, Rxn, Met, findempty,
                        av, va, nzabs_range, set_met!, set_rxn!,
                        struct_to_dict, isfixxed, ub, ub!, lb!, lb, FWD_SUFFIX, BKWD_SUFFIX

import Chemostat.SimulationUtils: tagprintln_inmw
import Chemostat.SteadyState: apply_bound!
import Chemostat.LP: fba, fva_preprocess


## -------------------------------------------------------------------
# Run add https://github.com/josePereiro/Chemostat_Heerden2013.jl in the Julia Pkg REPL to install the
# package, then you must activate the package enviroment (see README)
import Chemostat_Heerden2013: HeerdenData, BegData, iJR904, save_data, load_data
import Chemostat_Heerden2013.iJR904: OBJ_IDER, ATPM_IDER, COST_IDER
const Hd  = HeerdenData;
const Bd  = BegData
const iJR = iJR904


## -------------------------------------------------------------------
# BASE MODEL
# Load Mat file
println("Original .mat model")
src_file = iJR.MODEL_RAW_MAT_FILE
mat_model = MAT.matread(src_file)["model"]
model = MetNet(mat_model; reshape = true);
tagprintln_inmw("MAT MODEL LOADED", 
    "\nfile:             ", relpath(src_file), 
    "\nfile size:        ", filesize(src_file), " bytes", 
    "\nmodel size:       ", size(model),
    "\nnzabs_range:      ", nzabs_range(model.S),
)

## -------------------------------------------------------------------
# the maximal experimental growth rate in Heerden2013 is ~0.2 1/h
# The raw model present a growth rate bigger than that, so it is ok
# to use it directly as base model
fbaout = fba(model, OBJ_IDER);
tagprintln_inmw("FBA SOLUTION", 
    "\nobj_ider:         " , OBJ_IDER,
    "\nfba obj_val:      ", av(model, fbaout, OBJ_IDER),
    "\nmax exp obj_val:  ", maximum(Hd.val("D"))
)

## -------------------------------------------------------------------
# Set bounds
# The abs maximum bounds will be set to 100
abs_max_bound = 100
tagprintln_inmw("CLAMP BOUNDS", 
    "\nabs max bound: ", abs_max_bound
)
foreach(model.rxns) do ider
        isfixxed(model, ider) && return # fixxed reaction are untouched

        old_ub = ub(model, ider)
        new_ub = old_ub == 0.0 ? 0.0 : abs_max_bound
        ub!(model, ider, new_ub)

        old_lb = lb(model, ider)
        new_lb = old_lb == 0.0 ? 0.0 : -abs_max_bound
        lb!(model, ider, new_lb)
end

## -------------------------------------------------------------------
exchs = exchanges(model)
tagprintln_inmw("CLOSE EXCANGES", 
    "\nexchanges: ", exchs |> length
)
# Close, for now, all exchanges for avoiding it to be in revs
# The reversible reactions will be splited for modeling cost
# Exchanges have not associated cost, so, we do not split them
foreach(exchs) do idx
    ub!(model, idx, 0.0) # Closing all outtakes
    lb!(model, idx, 0.0) # Closing all intakes
end


## -------------------------------------------------------------------
# Enzymatic cost info
# The cost will be introduced as a reaction, we follow the same cost models as 
# Beg et al. (2007): https://doi.org/10.1073/pnas.0609845104.
# A new balance equations is then added:
#        Σ(rᵢ*costᵢ) + tot_cost = 0
#    Because the cost coefficients (costᵢ) < 0 (it resamble a reactant), the system must allocate 
#    the fluxes (rᵢ) so that Σ(rᵢ*costᵢ) = tot_cost, and tot_cost
#    are usually bounded [0.0, 1.0]
cost_info = Dict()
fwd_ider(rxn) = string(rxn, FWD_SUFFIX);
bkwd_ider(rxn) = string(rxn, BKWD_SUFFIX);
for rxn in model.rxns
    # The exchanges, the atpm and the biomass are synthetic reactions, so, 
    # they have should not have an associated enzimatic cost 
    any(startswith.(rxn, ["EX_", "DM_"])) && continue
    rxn == OBJ_IDER && continue
    rxn == ATPM_IDER && continue
        
    # Only the internal, non reversible reactions have an associated cost
    # We will split the rev reactions, so we save the cost for both versions (fwd, bkwd)
    if isrev(model, rxn)
        cost_info[fwd_ider(rxn)] = -iJR.beg_enz_cost(rxn)
        cost_info[bkwd_ider(rxn)] = -iJR.beg_enz_cost(rxn)
    else
        cost_info[rxn] = -iJR.beg_enz_cost(rxn)
    end
end

## -------------------------------------------------------------------
# Spliting revs
tagprintln_inmw("SPLITING REVS", 
    "\nfwd_suffix:      ", FWD_SUFFIX,
    "\nbkwd_suffix:     ", BKWD_SUFFIX,
)

model = split_revs(model;

    get_fwd_ider = fwd_ider,
    get_bkwd_ider = bkwd_ider,
);

## -------------------------------------------------------------------
# Adding cost raction
cost_met_id = "cost"
cost_exch_id = COST_IDER
tagprintln_inmw("ADDING COST", 
    "\ncosts to add: ", cost_info |> length,
    "\nmin abs coe:  ", cost_info |> values .|> abs |> minimum,
    "\nmax abs coe:  ", cost_info |> values .|> abs |> maximum,
    "\ncost met id:  ", cost_met_id,
    "\ncost exch id: ", cost_exch_id
)

M, N = size(model)
cost_met = Met(cost_met_id, S = collect(values(cost_info)), rxns = collect(keys(cost_info)), b = 0.0)
model = expanded_model(model, M + 1, N + 1)
set_met!(model, findempty(model, :mets), cost_met)
cost_exch = Rxn(cost_exch_id, S = [1.0], mets = [cost_met_id], lb = -abs_max_bound, ub = 0.0, c = 0.0)
set_rxn!(model, findempty(model, :rxns), cost_exch);

## -------------------------------------------------------------------
# Set base exchanges
tagprintln_inmw("SETTING EXCHANGES") 
# To control the intakes just the metabolites defined in the 
# base_intake_info (The minimum medium) will be opened.
# The base model will be constraint as in a cultivation with xi = 1.0
# see Cossios paper (see README)
ξ = 1.0 
# println("Minimum medium: ", iJR.base_intake_info)
foreach(exchs) do idx
    ub!(model, idx, abs_max_bound) # Opening all outakes
    lb!(model, idx, 0.0) # Closing all intakes
end

# see Cossios paper (see README) for details in the Chemostat bound constraint
apply_bound!(model, ξ, iJR.base_intake_info);

# tot_cost is the exchange that controls the bounds of the 
# enzimatic cost contraint, we bound it to [0, 1.0]
lb!(model, cost_exch_id, 0.0);
ub!(model, cost_exch_id, 1.0);

## -------------------------------------------------------------------
# Exch_met_map

# A quick way to get exchages from mets and te other way around
exch_met_map = Dict()
for exch_ in model.rxns[findall((id) -> any(startswith.(id, ["EX_", "DM_"])), model.rxns)]
    mets_ = model.mets[rxn_mets(model, exch_)]
    length(mets_) != 1 && continue
    exch_met_map[exch_] = mets_[1]
    exch_met_map[mets_[1]] = exch_
end;

# saving
save_data(iJR.EXCH_MET_MAP_FILE, exch_met_map)


## -------------------------------------------------------------------
## FVA PREPROCESSING
model = fva_preprocess(model, 
#     eps = 1-9, # This avoid blocking totally any reaction
    verbose = true);

fbaout = fba(model, OBJ_IDER, COST_IDER);
tagprintln_inmw("FBA SOLUTION", 
    "\nobj_ider:         " , OBJ_IDER,
    "\nfba obj_val:      ", av(model, fbaout, OBJ_IDER),
    "\nmax exp obj_val:  ", maximum(Hd.val("D")),
    "\ncost_ider:         " , COST_IDER,
    "\nfba cost_val:      ", av(model, fbaout, COST_IDER),
)
Chemostat.Utils.summary(model, fbaout)


## -------------------------------------------------------------------
# Saving model
dict_model = model |> struct_to_dict |> compress_dict
save_data(iJR.BASE_MODEL_FILE, dict_model);


