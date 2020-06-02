module iJR904

import ..Chemostat_Heerden2013: PROJ_ROOT, RAW_DATA_DIR, PROCESSED_DATA_DIR
# import Chemostat.Utils: MetNet
import ..BegData

include("dirs_and_files.jl")
include("maps.jl")
include("enz_cost.jl")
include("base_intake_info.jl")

end