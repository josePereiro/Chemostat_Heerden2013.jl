# Beg et al. (2007): https://doi.org/10.1073/pnas.0609845104.
module BegData
    import ..Chemostat_Heerden2013: PROJ_ROOT, RAW_DATA_DIR, PROCESSED_DATA_DIR
    import CSV

    include("dirs_and_files.jl")
    include("enz_data.jl")
    
end