module iJR904

import BSON
import ..Chemostat_Heerden2013: PROJ_ROOT, DATA_DIR, FIGURES_DATA_DIR, RAW_DATA_DIR, PROCESSED_DATA_DIR
import ..Chemostat_Heerden2013.Chemostat.Utils: load_data
import ..BegData
import CSV
import DataFrames: DataFrame

include("const.jl")
include("dirs_and_files.jl")
include("maps.jl")
include("beg_enz_cost.jl")
include("base_intake_info.jl")

function __init__()
    _create_dirs()
end

end