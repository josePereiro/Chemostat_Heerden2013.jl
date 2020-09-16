module iJR904

import BSON
import ..Chemostat_Heerden2013: PROJ_ROOT, DATA_DIR, FIGURES_DATA_DIR, RAW_DATA_DIR, PROCESSED_DATA_DIR, load_data
import ..BegData
import CSV
import DataFrames: DataFrame

include("dirs_and_files.jl")
include("maps.jl")
include("beg_enz_cost.jl")
include("base_intake_info.jl")
include("const.jl")

function __init__()
    load_exch_met_map()
end

end