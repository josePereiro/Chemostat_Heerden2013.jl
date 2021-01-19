import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

import CSV
import MAT
using DataFrames
using Serialization

## -------------------------------------------------------------------
# Run add https://github.com/josePereiro/Chemostat_Heerden2013.jl in the Julia Pkg REPL to install the
# package, then you must activate the package enviroment (see README)
import Chemostat_Heerden2013
const Hd  = Chemostat_Heerden2013.HeerdenData;
const ECC = Chemostat_Heerden2013.EColiCore

## -------------------------------------------------------------------
# run add "https://github.com/josePereiro/Chemostat" in the 
# julia Pkg REPL for installing the package
import Chemostat
const ChU = Chemostat.Utils
const ChSS = Chemostat.SteadyState
const ChLP = Chemostat.LP

## -------------------------------------------------------------------
# Tools
function partial_test(model, title  = "PARTIAL TEST")
    fbaout = ChLP.fba(model, ECC.BIOMASS_IDER);
    ChU.tagprintln_inmw(title, 
        "\nsize:             ", size(model),
        "\nobj_ider:         " , ECC.BIOMASS_IDER,
        "\nfba obj_val:      ", ChU.av(model, fbaout, ECC.BIOMASS_IDER),
        "\nmax exp obj_val:  ", maximum(Hd.val("D")),
    )
end

## -------------------------------------------------------------------
# BASE MODEL
# Load Mat file
println("Original .mat model")
src_file = ECC.MODEL_RAW_MAT_FILE
mat_model = MAT.matread(src_file)["e_coli_core"]
model = ChU.MetNet(mat_model; reshape = true);

## -------------------------------------------------------------------
# the maximal experimental growth rate in Heerden2013 is ~0.2 1/h
# The raw model present a growth rate bigger than that, so it is ok
# to use it directly as base model
partial_test(model)

## -------------------------------------------------------------------
# Set bounds
# The abs maximum bounds will be set to ECC.MAX_ABS_BOUND
ChU.tagprintln_inmw("CLAMP BOUNDS", 
    "\nabs max bound: ", ECC.MAX_ABS_BOUND
)
ChU.clampfields!(model, [:lb, :ub]; 
    abs_max = ECC.MAX_ABS_BOUND, zeroth = 0.0)
partial_test(model)

## -------------------------------------------------------------------
# Exchanges
exchs = ChU.exchanges(model)
ChU.tagprintln_inmw("SETTING EXCHANGES") 
# To control the intakes just the metabolites defined in the 
# base_intake_info (The minimum medium) will be opened.
# The base model will be constraint as in a cultivation with xi = 1.0
# see Cossios paper (see README)
ξ = 1.0 
# println("Minimum medium: ", ECC.base_intake_info)
foreach(exchs) do idx
    ChU.ub!(model, idx, ECC.MAX_ABS_BOUND) # Opening all outakes
    ChU.lb!(model, idx, 0.0) # Closing all intakes
end

# see Cossios paper (see README) for details in the Chemostat bound constraint
ChSS.apply_bound!(model, ξ, ECC.base_intake_info);

## -------------------------------------------------------------------
# Eliminate similars
ChU.tagprintln_inmw("ELIMINATE SIMILAR RXNS") 
ChU.summary.([model], ["SUCDi", "FRD7"])
ChU.bounds!(model, "SUCDi", 0.0, 0.0)
ChU.bounds!(model, "FRD7", -ECC.MAX_ABS_BOUND, ECC.MAX_ABS_BOUND)
ChU.summary.([model], ["SUCDi", "FRD7"])

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
ChU.save_data(ECC.EXCH_MET_MAP_FILE, exch_met_map)

## -------------------------------------------------------------------
partial_test(model)

# FVA PREPROCESSING
model = ChLP.fva_preprocess(model, 
    batchlen = 1,
    check_obj = ECC.BIOMASS_IDER,
    verbose = true
)
partial_test(model)

## -------------------------------------------------------------------
ChU.summary(model)

## -------------------------------------------------------------------
# Saving model
dict_model = model |> ChU.struct_to_dict
ChU.save_data(ECC.BASE_MODELS_FILE, dict_model);