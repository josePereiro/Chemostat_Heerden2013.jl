import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

@time begin
    import CSV
    import MAT
    using DataFrames
    import SparseArrays
    using Serialization

    # -------------------------------------------------------------------
    # run add "https://github.com/josePereiro/Chemostat" in the 
    # julia Pkg REPL for installing the package
    import Chemostat
    const ChU = Chemostat.Utils
    const ChSS = Chemostat.SteadyState
    const ChLP = Chemostat.LP

    # -------------------------------------------------------------------
    # Run add https://github.com/josePereiro/Chemostat_Heerden2013.jl in the Julia Pkg REPL to install the
    # package, then you must activate the package enviroment (see README)
    import Chemostat_Heerden2013
    const ChH = Chemostat_Heerden2013
    const Hd  = ChH.HeerdenData;
    const Bd  = ChH.BegData
    const iJR = ChH.iJR904
end

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
src_file = iJR.rawdir("iJR904.mat")
mat_model = MAT.matread(src_file)["model"]
model = ChU.MetNet(mat_model; reshape = true);

## -------------------------------------------------------------------
# the maximal experimental growth rate in Heerden2013 is ~0.2 1/h
# The raw model present a growth rate bigger than that, so it is ok
# to use it directly as base model
partial_test(model)

## -------------------------------------------------------------------
# Gene intactivations
# Heerden2013 describes some gene modifications
# https://doi.org/10.1186/1475-2859-12-80. Table 2
ChU.tagprintln_inmw("GENE MODIFICATIONS")
for (name, iders) in iJR.load_Hd_to_inactivate_map()

    name in ["Aspartate aminotransferase"] && continue ## Closing this one are lethal

    ChU.bounds!.([model], iders, 0.0, 0.0)
    @info("Closing gene modificated", iders)
    partial_test(model)
    println()
end

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
let 
    ξ = 1.0 
    base_intake_info = iJR.load_base_intake_info()

    # println("Minimum medium: ", iJR.base_intake_info)
    foreach(exchs) do idx
        ChU.ub!(model, idx, iJR.MAX_ABS_BOUND) # Opening all outakes
        ChU.lb!(model, idx, 0.0) # Closing all intakes
    end

    # see Cossios paper (see README) for details in the Chemostat bound constraint
    ChSS.apply_bound!(model, ξ, base_intake_info);
    # tot_cost is the exchange that controls the bounds of the 
    # enzimatic cost contraint, we bound it to [0, 1.0]
    ChU.lb!(model, cost_exch_id, 0.0);
    ChU.ub!(model, cost_exch_id, 1.0);
end

## -------------------------------------------------------------------
function scale_model(model, scale_factor)
    base_nzabs_range = ChU.nzabs_range(model.S)
    base_size = size(model)
    
    # Scale model (reduce S ill-condition)
    model = ChU.well_scaled_model(model, scale_factor)
    
    scl_size = size(model)
    scl_nzabs_range = ChU.nzabs_range(model.S)

    @info("Model", exp, 
        base_size, base_nzabs_range, 
        scl_size, scl_nzabs_range
    ); println()
    return model
end

## -------------------------------------------------------------------
# FVA PREPROCESSING
compressed(model) = model |> ChU.struct_to_dict |> ChU.compressed_copy
MODEL_FILE = iJR.procdir("base_models.bson")
BASE_MODELS = isfile(MODEL_FILE) ? 
    ChU.load_data(MODEL_FILE) : 
    Dict("load_model" => compressed(model))
cGLCs = Hd.val("cGLC")
for (exp, cGLC) in enumerate(cGLCs)

    D = get!(BASE_MODELS, "fva_models", Dict())
    ChU.tagprintln_inmw("DOING FVA", 
        "\nexp:             ", exp,
        "\ncGLC:            ", cGLC,
        "\ncProgress:       ", length(D),
        "\n"
    )
    haskey(D, exp) && continue # cached

    ## -------------------------------------------------------------------
    # prepare model
    # Scale model (reduce S ill-condition)
    scale_factor = 1000.0
    model0 = scale_model(deepcopy(model), scale_factor)
    intake_info = iJR.load_base_intake_info()
    # The only difference between experiments is the feed medium 
    # concentration.
    ξ = Hd.val("xi", exp)
    intake_info[iJR.EX_GLC_IDER]["c"] = cGLC
    ChSS.apply_bound!(model0, ξ, intake_info; 
        emptyfirst = true, ignore_miss = true
    )
        
    ## -------------------------------------------------------------------
    # fva
    partial_test(model0)
    fva_model = ChLP.fva_preprocess(model0, 
        batchlen = 50,
        check_obj = iJR.BIOMASS_IDER,
        verbose = true
    )
    partial_test(fva_model)

    # storing
    D[exp] = compressed(fva_model)
end

## -------------------------------------------------------------------
# MAX MODEL
let
    # This model is bounded by the maximum rates found for EColi.
    # Data From:
    # Varma, (1993): 2465–73. https://doi.org/10.1128/AEM.59.8.2465-2473.1993.
    # Extract max exchages from FIG 3 to form the maximum polytope

    ChU.tagprintln_inmw("DOING MAX MODEL", 
        "\n"
    )
    
    scale_factor = 1000.0
    max_model = scale_model(deepcopy(model), scale_factor)
    
    # Biomass
    # 2.2 1/ h
    ChU.bounds!(max_model, iJR.BIOMASS_IDER, 0.0, 2.2)
    
    Hd_rxns_map = iJR.load_rxns_map() 
    # 40 mmol / gDW h
    ChU.bounds!(max_model, Hd_rxns_map["GLC"], -40.0, 0.0)
    # 45 mmol/ gDW
    ChU.bounds!(max_model, Hd_rxns_map["AC"], 0.0, 40.0)
    # 55 mmol/ gDW h
    ChU.bounds!(max_model, Hd_rxns_map["FORM"], 0.0, 55.0)
    # 20 mmol/ gDW h
    ChU.bounds!(max_model, Hd_rxns_map["O2"], -20.0, 0.0)
    
    # fva
    max_model = ChLP.fva_preprocess(max_model, 
        check_obj = iJR.BIOMASS_IDER,
        verbose = true
    );

    ## -------------------------------------------------------------------
    test_model = deepcopy(max_model)
    for (exp, D) in Hd.val(:D) |> enumerate
        cgD_X = Hd.cval(:GLC, exp) * Hd.val(:D, exp) / Hd.val(:DCW, exp)
        ChU.lb!(test_model, iJR.EX_GLC_IDER, -cgD_X)
        fbaout = ChLP.fba(test_model, iJR.BIOMASS_IDER, iJR.COST_IDER)
        biom = ChU.av(test_model, fbaout, iJR.BIOMASS_IDER)
        cost = ChU.av(test_model, fbaout, iJR.COST_IDER)
        @info("Test", exp, cgD_X, D, biom, cost); println()
    end
    
    ## -------------------------------------------------------------------
    # saving
    BASE_MODELS["max_model"] = compressed(max_model)
end;

## -------------------------------------------------------------------
# SAVING
ChU.save_data(MODEL_FILE, BASE_MODELS);
