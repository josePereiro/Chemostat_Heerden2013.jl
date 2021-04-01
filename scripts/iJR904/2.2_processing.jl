import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid

    # -------------------------------------------------------------------
    # Run add https://github.com/josePereiro/Chemostat_Heerden2013.jl in 
    # the Julia Pkg REPL to install the package, then you must activate 
    # the package enviroment (see README)
    import Chemostat_Heerden2013
    const ChHd = Chemostat_Heerden2013
    const Hd  = ChHd.HeerdenData
    const BD  = ChHd.BegData
    const iJR = ChHd.iJR904

    # -------------------------------------------------------------------
    # run add "https://github.com/josePereiro/Chemostat" in the 
    # julia Pkg REPL for installing the package
    import Chemostat
    import Chemostat.LP.MathProgBase

    const Ch = Chemostat
    
    const ChU = Ch.Utils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP
    const ChEP = Ch.MaxEntEP
    const ChSU = Ch.SimulationUtils

    import ChemostatPlots
    const ChP = ChemostatPlots

    import UtilsJL
    const UJL = UtilsJL

    using Statistics
    using ProgressMeter
    using Base.Threads
    using Serialization

    # -------------------------------------------------------------------
    using Plots, FileIO
    import GR
    GR.inline("png")
end

## ----------------------------------------------------------------------------
# ME data index
ME_INDEX = ChU.load_data(iJR.MAXENT_VARIANTS_ME_INDEX_FILE; verbose = false);

## ----------------------------------------------------------------------------
# LP data
LPDAT = ChU.load_data(iJR.LP_DAT_FILE; verbose = false)

# -------------------------------------------------------------------
const ME_Z_OPEN_G_OPEN          = :ME_Z_OPEN_G_OPEN           # Do not use extra constraints
const ME_MAX_POL                = :ME_MAX_POL                 # Match ME and Dy biom average and constraint av_ug
const ME_Z_EXPECTED_G_BOUNDED   = :ME_Z_EXPECTED_G_BOUNDED    # Match ME and Dy biom average and constraint av_ug
const ME_Z_EXPECTED_G_MOVING    = :ME_Z_EXPECTED_G_MOVING     # 
const ME_Z_FIXXED_G_BOUNDED     = :ME_Z_FIXXED_G_BOUNDED      # Fix biom around observed

# -------------------------------------------------------------------
fileid = "2.1"
function mysavefig(p, pname; params...) 
    fname = UJL.mysavefig(p, string(fileid, "_", pname), iJR.MODEL_FIGURES_DIR; params...)
    @info "Plotting" fname
end
myminmax(a::Vector) = (minimum(a), maximum(a))
CONC_IDERS = ["GLC", "SUCC", "AC", "FORM"]
FLX_IDERS = ["GLC", "SUCC", "AC", "FORM"]

ALL_METHODS = [
    # ME_Z_OPEN_G_OPEN, 
    ME_MAX_POL, 
    # ME_Z_FIXXED_G_BOUNDED, 
    # ME_Z_EXPECTED_G_BOUNDED, 
    # ME_Z_EXPECTED_G_MOVING
]

Hd_rxns_map = iJR.load_Hd_rxns_map()
Hd_mets_map = iJR.load_Hd_mets_map()

## ----------------------------------------------------------------------------
# collect DAT
include("2.2.1_collect_DAT.jl")

## ----------------------------------------------------------------------------
EXPS = DAT[:EXPS]

exp_colors = let
    colors = Plots.distinguishable_colors(length(EXPS))
    Dict(exp => color for (exp, color) in zip(EXPS, colors))
end

ider_colors = Dict(
    "GLC" => :red, "SUCC" => :yellow,
    "AC" => :orange, "FORM" => :blue,
    "D" => :black,
)

method_colors = Dict(
    ME_Z_OPEN_G_OPEN => :red,
    ME_MAX_POL => :blue,
    ME_Z_EXPECTED_G_BOUNDED => :orange,
    ME_Z_EXPECTED_G_MOVING => :purple,
    ME_Z_FIXXED_G_BOUNDED => :green,
)

## -------------------------------------------------------------------
# PLOTS
# -------------------------------------------------------------------

## -------------------------------------------------------------------
# polytope box volume
include("2.2.2_pol_vox_volume.jl")

## -------------------------------------------------------------------
# Var study
include("2.2.3_var_study.jl")

## -------------------------------------------------------------------
# MSE study
include("2.2.4_MSE_study.jl")

## -------------------------------------------------------------------
# proj 2D
include("2.2.5_proj2d_study.jl")

## -------------------------------------------------------------------
# Beta tendencies
include("2.2.6_Beta_tendencies.jl")

## -------------------------------------------------------------------
# flx correlations
include("2.2.7_flx_correlations.jl")

## -------------------------------------------------------------------
# inbounds study
include("2.2.7_flx_correlations.jl")

## -------------------------------------------------------------------
# marginals
include("2.2.9_marginals.jl")

## -------------------------------------------------------------------
# legends
include("2.2.10_leyends.jl")