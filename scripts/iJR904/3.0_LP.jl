import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid, SpinLock

    #  ----------------------------------------------------------------------------
    # Run add https://github.com/josePereiro/Chemostat_Heerden2013.jl in the Julia Pkg REPL to install the
    # package, then you must activate the package enviroment (see README)
    import Chemostat_Heerden2013
    const ChHd = Chemostat_Heerden2013
    const Hd  = ChHd.HeerdenData;
    const BD  = ChHd.BegData;
    const iJR = ChHd.iJR904

    #  ----------------------------------------------------------------------------
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

    import JuMP, GLPK
    using Serialization
    import UtilsJL
    const UJL = UtilsJL
end

## -----------------------------------------------------------------------------------------------
const FBA_MIN_COST = :FBA_MIN_COST
const FBA_OPEN = :FBA_OPEN
const YIELD = :YIELD
const FBA_MAX_VG_YIELD = :FBA_MAX_VG_YIELD
const EXPS = 1:13

## -----------------------------------------------------------------------------------------------
function base_model(exp)
    BASE_MODELS = ChU.load_data(iJR.BASE_MODELS_FILE; verbose = false);
    model_dict = BASE_MODELS["fva_models"][exp]
    model = ChU.MetNet(;model_dict...) |> ChU.uncompressed_model
end

## -----------------------------------------------------------------------------------------------
# Data container
LPDAT = UJL.DictTree()

## -----------------------------------------------------------------------------------------------
# YIELD
# include("3.0.1_YIELD.jl")

## -------------------------------------------------------------------
# FBA
include("3.0.2_FBA.jl")

## -------------------------------------------------------------------
ChU.save_data(iJR.LP_DAT_FILE, LPDAT)
