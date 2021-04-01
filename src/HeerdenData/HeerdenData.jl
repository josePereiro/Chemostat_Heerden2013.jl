# Heerden2013 https=>//doi.org/10.1186/1475-2859-12-80.

module HeerdenData
    import ..Chemostat_Heerden2013
    const ChH = Chemostat_Heerden2013
    import CSV
    import DataFrames: DataFrame

    import UtilsJL
    const UJL = UtilsJL
    UJL.gen_sub_proj(@__MODULE__)

    include("data_interface.jl")

    function __init__()
        _load_cont_cul_data()
        UJL.create_proj_dirs(@__MODULE__)
    end

end