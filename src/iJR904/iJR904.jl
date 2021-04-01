module iJR904

    import BSON
    import ..Chemostat_Heerden2013
    const ChH = Chemostat_Heerden2013
    import Chemostat
    const ChU = Chemostat.Utils
    import ..BegData
    const Bd = BegData
    import CSV
    import DataFrames: DataFrame

    import UtilsJL
    const UJL = UtilsJL
    UJL.gen_sub_proj(@__MODULE__)

    include("const.jl")
    include("maps.jl")
    include("beg_enz_cost.jl")
    include("base_intake_info.jl")
    include("load_model.jl")

    function __init__()
        UJL.create_proj_dirs(@__MODULE__)
    end

end