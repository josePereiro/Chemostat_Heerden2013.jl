# Here I include some maps between the experimental
# data ids and the model ids

## -------------------------------------------------------------------
# maps between Heerden2013 https://doi.org/10.1186/1475-2859-12-80 
# data and the model

# Mets map
const Hd_mets_map = Dict();
function _set_Hd_mets_map()
    empty!(Hd_mets_map)
    Hd_mets_map["GLC"] = "glc__D_e";
    Hd_mets_map["SA"] = "succ_e";
    Hd_mets_map["AcA"] = "ac_e";
    Hd_mets_map["FA"] = "for_e";
    Hd_mets_map["MA"] = "mal__L_e";
    for (k, v) in Hd_mets_map
        Hd_mets_map[v] = k
    end
    Hd_mets_map
end

# base model exch met map
# A quick way to get exchages from mets and te other way around
const exch_met_map = Dict()
function _load_exch_met_map()
    !isfile(EXCH_MET_MAP_FILE) && return exch_met_map
    empty!(exch_met_map)
    merge!(exch_met_map, load_data(EXCH_MET_MAP_FILE; verbose = false))
    return exch_met_map
end
