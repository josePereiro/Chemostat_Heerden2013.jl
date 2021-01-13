import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

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
ChU.set_cache_dir(iJR.CACHE_DIR)

## ----------------------------------------------------------------------------
let
    model = ChU.MetNet(ChU.load_data(iJR.BASE_MODEL_FILE))
    model = ChU.uncompressed_model(model)
    objidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
    M, N = size(model)
end