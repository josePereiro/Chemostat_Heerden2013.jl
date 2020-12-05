module EColiCore

import ..Chemostat_Heerden2013: PROJ_ROOT, DATA_DIR, FIGURES_DATA_DIR, RAW_DATA_DIR, PROCESSED_DATA_DIR
import ..Chemostat_Heerden2013.Chemostat.Utils: load_data

include("const.jl")
include("file_and_dirs.jl")
include("base_intake_info.jl")
include("maps.jl")

function __init__()
    _create_dirs()
    _set_Hd_mets_map()
    _load_exch_met_map()
end

end