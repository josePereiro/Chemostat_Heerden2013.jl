# (GEM) iJR904 (Reed et al., 2003) (download link: https://darwin.di.uminho.pt/models)
# TODO: change to alternative link!!!
# (GEM) Alternative link http://bigg.ucsd.edu/static/models/iJR904.mat
module Chemostat_Heerden2013

    import Chemostat
    const Ch = Chemostat

    import UtilsJL
    const UJL = UtilsJL
    UJL.gen_top_proj(@__MODULE__)

    include("Utils/Utils.jl")
    include("BegData/BegData.jl")
    include("HeerdenData/HeerdenData.jl")
    include("iJR904/iJR904.jl")
    # include("EColiCore/EColiCore.jl")

    function __init__()
        UJL.create_proj_dirs(@__MODULE__)
    end

end # module
