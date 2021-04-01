# # DIRS
# const HEERDEN_RAW_DATA_DIR = 
#     joinpath(RAW_DATA_DIR, "Heerden2013")
# const HEERDEN_PROCESSED_DATA_DIR = 
#     joinpath(PROCESSED_DATA_DIR, basename(HEERDEN_RAW_DATA_DIR))
# const HEERDEN_FIGURES_DIR = 
#     joinpath(FIGURES_DATA_DIR, basename(HEERDEN_RAW_DATA_DIR))

# function _create_dirs()
#     for dir in [HEERDEN_PROCESSED_DATA_DIR, HEERDEN_FIGURES_DIR]
#         try; mkpath(dir); catch err; @warn("Error creating dir", dir, err); end
#     end
# end
     

# # FILES
# const HEERDEN_CONT_CUL_DATA_ORIG_FILE = 
#     joinpath(HEERDEN_RAW_DATA_DIR, "heerden2013___cont_cult_data.tsv")
# const HEERDEN_CONT_CUL_DATA_CONV_FILE = 
#     joinpath(HEERDEN_PROCESSED_DATA_DIR, basename(HEERDEN_CONT_CUL_DATA_ORIG_FILE))
# const HEERDEN_GENE_MODS_FILE = 
#     joinpath(HEERDEN_RAW_DATA_DIR, "heerden2013___data", "heerden2013___gene_modifications.tsv")
