# I map some model ids to the enzimatic cost data reported in 
# Beg2007, so, I will set the cost of this reactions using this data
# The rest of the model reactions will get the average enzimatic cost also 
# reported by Beg.
_beg_enz_cost_dict = Dict()
for (i, beg_rxn) in enumerate(BegData.beg_enz_data.Enzyme)
    for model_rnx in beg_rxns_map[beg_rxn]
        if model_rnx != ""
            _beg_enz_cost_dict[model_rnx] = BegData.beg_enz_data[i, Symbol("ai/xi (h g / mmol)")]
        end
    end
end

function beg_enz_cost(model_rxn)
    return get(_beg_enz_cost_dict, model_rxn, BegData.ave_a)
end