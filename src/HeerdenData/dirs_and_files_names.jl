# DIRS
const HEERDEN_RAW_DATA_DIR = joinpath(RAW_DATA_DIR, "heerden2013___data")
const HEERDEN_PROCESSED_DATA_DIR = joinpath(PROCESSED_DATA_DIR, basename(HEERDEN_RAW_DATA_DIR))
mkpath(HEERDEN_PROCESSED_DATA_DIR)    

# FILES
const HEERDEN_CONT_CUL_DATA_ORIG_FILE = joinpath(HEERDEN_RAW_DATA_DIR, "heerden2013___cont_cult_data.tsv")
const HEERDEN_CONT_CUL_DATA_CONV_FILE = joinpath(HEERDEN_PROCESSED_DATA_DIR, basename(HEERDEN_CONT_CUL_DATA_ORIG_FILE))
Base.include_dependency(HEERDEN_CONT_CUL_DATA_CONV_FILE)
const HEERDEN_GENE_MODS_FILE = joinpath(HEERDEN_RAW_DATA_DIR, "heerden2013___data/heerden2013___gene_modifications.tsv")
