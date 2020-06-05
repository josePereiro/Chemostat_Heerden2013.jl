# DIRS
const MODEL_RAW_DIR = joinpath(RAW_DATA_DIR, "iJR904")
const MODEL_PROCESSED_DATA_DIR = joinpath(PROCESSED_DATA_DIR, basename(MODEL_RAW_DIR))
mkpath(MODEL_PROCESSED_DATA_DIR)
const MODEL_FIGURES_DIR = joinpath(FIGURES_DATA_DIR, basename(MODEL_RAW_DIR))
mkpath(MODEL_FIGURES_DIR)

# FILES
const MODEL_RAW_MAT_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "iJR904.mat")
Base.include_dependency(MODEL_RAW_MAT_FILE)
const BASE_MODEL_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "iJR904___base_model.jls")
Base.include_dependency(BASE_MODEL_FILE)
const EXCH_MET_MAP_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "exch_met_map.csv")
Base.include_dependency(EXCH_MET_MAP_FILE)