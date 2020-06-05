module iJR904

import ..Chemostat_Heerden2013: PROJ_ROOT, DATA_DIR, FIGURES_DATA_DIR, RAW_DATA_DIR, PROCESSED_DATA_DIR
# import Chemostat.Utils: MetNet
import ..BegData
import CSV
import DataFrames: DataFrame

include("dirs_and_files.jl")
include("maps.jl")
include("beg_enz_cost.jl")
include("base_intake_info.jl")

end