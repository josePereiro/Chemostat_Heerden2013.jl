# Heerden2013 https=>//doi.org/10.1186/1475-2859-12-80.

module HeerdenData
    import ..Chemostat_Heerden2013: PROJ_ROOT, RAW_DATA_DIR, PROCESSED_DATA_DIR
    import CSV
    import DataFrames: DataFrame

    include("dirs_and_files_names.jl")
    include("data_interface.jl")

    function __init__()
        _load_cont_cul_data()
    end

end