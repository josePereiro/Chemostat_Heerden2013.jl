# DIRS
const MODEL_RAW_DIR = joinpath(RAW_DATA_DIR, PROJ_IDER)
const MODEL_PROCESSED_DATA_DIR = joinpath(PROCESSED_DATA_DIR, basename(MODEL_RAW_DIR))
const MODEL_FIGURES_DIR = joinpath(FIGURES_DATA_DIR, basename(MODEL_RAW_DIR))
const CACHE_DIR = joinpath(MODEL_PROCESSED_DATA_DIR, "cache")

function _create_dirs()
    for dir in [MODEL_PROCESSED_DATA_DIR, 
                MODEL_FIGURES_DIR, 
                CACHE_DIR
            ]
        isdir(dir) && continue
        try; mkpath(dir); catch end
    end
end

# FILES
const MODEL_RAW_MAT_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "iJR904.mat")
const BASE_MODELS_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "base_models.bson")
const EXCH_MET_MAP_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "exch_met_map.bson")
const MAXENT_B0SEEDS_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "maxent_b0seed_bundles.bson")
const MAXENT_VARIANTS_INDEX_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "maxent_variants_index.bson")
const CORR_DAT_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "corr_dat_file.bson")
const LP_DAT_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "lp_dat.bson")