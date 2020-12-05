import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

import CSV
import MAT
using DataFrames
using Serialization

## -------------------------------------------------------------------
# run add "https://github.com/josePereiro/Chemostat" in the 
# julia Pkg REPL for installing the package
import Chemostat
const ChU = Chemostat.Utils
const ChSS = Chemostat.SteadyState
const ChLP = Chemostat.LP

## -------------------------------------------------------------------
# Run add https://github.com/josePereiro/Chemostat_Heerden2013.jl in the Julia Pkg REPL to install the
# package, then you must activate the package enviroment (see README)
import Chemostat_Heerden2013: HeerdenData, BegData, iJR904
const Hd  = HeerdenData;
const Bd  = BegData
const iJR = iJR904


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

## -------------------------------------------------------------------
# BASE MODEL
# Load Mat file
println("Original .mat model")
src_file = iJR.MODEL_RAW_MAT_FILE
mat_model = MAT.matread(src_file)["model"]
model = ChU.MetNet(mat_model; reshape = true);

## -------------------------------------------------------------------
# the maximal experimental growth rate in Heerden2013 is ~0.2 1/h
# The raw model present a growth rate bigger than that, so it is ok
# to use it directly as base model
partial_test(model)

## -------------------------------------------------------------------
# Set bounds
# The abs maximum bounds will be set to iJR.MAX_ABS_BOUND
ChU.tagprintln_inmw("CLAMP BOUNDS", 
    "\nabs max bound: ", iJR.MAX_ABS_BOUND
)
foreach(model.rxns) do ider
        ChU.isfixxed(model, ider) && return # fixxed reaction are untouched

        old_ub = ChU.ub(model, ider)
        new_ub = old_ub == 0.0 ? 0.0 : iJR.MAX_ABS_BOUND
        ChU.ub!(model, ider, new_ub)

        old_lb = ChU.lb(model, ider)
        new_lb = old_lb == 0.0 ? 0.0 : -iJR.MAX_ABS_BOUND
        ChU.lb!(model, ider, new_lb)
end
partial_test(model)

## -------------------------------------------------------------------
exchs = ChU.exchanges(model)
ChU.tagprintln_inmw("CLOSE EXCANGES", 
    "\nChU.exchanges: ", exchs |> length
)
# Close, for now, all ChU.exchanges for avoiding it to be in revs
# The reversible reactions will be splited for modeling cost
# Exchanges have not associated cost, so, we do not split them
foreach(exchs) do idx
    ChU.ub!(model, idx, 0.0) # Closing all outtakes
    ChU.lb!(model, idx, 0.0) # Closing all intakes
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
fwd_ider(rxn) = string(rxn, ChU.FWD_SUFFIX);
bkwd_ider(rxn) = string(rxn, ChU.BKWD_SUFFIX);
for rxn in model.rxns
    # The ChU.exchanges, the atpm and the biomass are synthetic reactions, so, 
    # they have should not have an associated enzimatic cost 
    any(startswith.(rxn, ["EX_", "DM_"])) && continue
    rxn == iJR.BIOMASS_IDER && continue
    rxn == iJR.ATPM_IDER && continue
        
    # Only the internal, non reversible reactions have an associated cost
    # We will split the rev reactions, so we save the cost for both versions (fwd, bkwd)
    if ChU.isrev(model, rxn)
        cost_info[fwd_ider(rxn)] = -iJR.beg_enz_cost(rxn)
        cost_info[bkwd_ider(rxn)] = -iJR.beg_enz_cost(rxn)
    else
        cost_info[rxn] = -iJR.beg_enz_cost(rxn)
    end
end

## -------------------------------------------------------------------
# Spliting revs
ChU.tagprintln_inmw("SPLITING REVS", 
    "\nfwd_suffix:      ", ChU.FWD_SUFFIX,
    "\nbkwd_suffix:     ", ChU.BKWD_SUFFIX,
)

model = ChU.split_revs(model;

    get_fwd_ider = fwd_ider,
    get_bkwd_ider = bkwd_ider,
);

## -------------------------------------------------------------------
# Adding cost raction
cost_met_id = "cost"
cost_exch_id = iJR.COST_IDER
ChU.tagprintln_inmw("ADDING COST", 
    "\ncosts to add: ", cost_info |> length,
    "\nmin abs coe:  ", cost_info |> values .|> abs |> minimum,
    "\nmax abs coe:  ", cost_info |> values .|> abs |> maximum,
    "\ncost met id:  ", cost_met_id,
    "\ncost exch id: ", cost_exch_id
)

M, N = size(model)
cost_met = ChU.Met(cost_met_id, S = collect(values(cost_info)), 
    rxns = collect(keys(cost_info)), b = 0.0)
model = ChU.expanded_model(model, M + 1, N + 1)
ChU.set_met!(model, ChU.findempty(model, :mets), cost_met)
cost_exch = ChU.Rxn(cost_exch_id, S = [1.0], mets = [cost_met_id], 
    lb = -iJR.MAX_ABS_BOUND, ub = 0.0, c = 0.0)
ChU.set_rxn!(model, ChU.findempty(model, :rxns), cost_exch);

## -------------------------------------------------------------------
# Set base ChU.exchanges
ChU.tagprintln_inmw("SETTING EXCHANGES") 
# To control the intakes just the metabolites defined in the 
# base_intake_info (The minimum medium) will be opened.
# The base model will be constraint as in a cultivation with xi = 1.0
# see Cossios paper (see README)
ξ = 1.0 
# println("Minimum medium: ", iJR.base_intake_info)
foreach(exchs) do idx
    ChU.ub!(model, idx, iJR.MAX_ABS_BOUND) # Opening all outakes
    ChU.lb!(model, idx, 0.0) # Closing all intakes
end

# see Cossios paper (see README) for details in the Chemostat bound constraint
ChSS.apply_bound!(model, ξ, iJR.base_intake_info);

# tot_cost is the exchange that controls the bounds of the 
# enzimatic cost contraint, we bound it to [0, 1.0]
ChU.lb!(model, cost_exch_id, 0.0);
ChU.ub!(model, cost_exch_id, 1.0);

## -------------------------------------------------------------------
# Exch_met_map

# A quick way to get exchages from mets and te other way around
exch_met_map = Dict()
for exch_ in model.rxns[findall((id) -> any(startswith.(id, ["EX_", "DM_"])), model.rxns)]
    mets_ = model.mets[ChU.rxn_mets(model, exch_)]
    length(mets_) != 1 && continue
    exch_met_map[exch_] = mets_[1]
    exch_met_map[mets_[1]] = exch_
end;

# saving
ChU.save_data(iJR.EXCH_MET_MAP_FILE, exch_met_map)

## -------------------------------------------------------------------
partial_test(model)

# FVA PREPROCESSING
model = ChLP.fva_preprocess(model, 
    batchlen = 50,
    check_obj = iJR.BIOMASS_IDER,
    verbose = true
)
partial_test(model)

## -------------------------------------------------------------------
ChU.summary(model)

## -------------------------------------------------------------------
# Saving model
dict_model = model |> ChU.struct_to_dict |> ChU.compressed_copy
ChU.save_data(iJR.BASE_MODEL_FILE, dict_model);