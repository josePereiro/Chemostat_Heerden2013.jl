import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

@time begin
    import CSV
    import MAT
    using DataFrames
    import SparseArrays
    using Serialization
    using ProgressMeter
    using UtilsJL
    const UJL = UtilsJL
    using Plots
    using Base.Threads

    # -------------------------------------------------------------------
    # run add "https://github.com/josePereiro/Chemostat" in the 
    # julia Pkg REPL for installing the package
    import Chemostat
    const ChU = Chemostat.Utils
    const ChSS = Chemostat.SteadyState
    const ChLP = Chemostat.LP

    # -------------------------------------------------------------------
    # Run add https://github.com/josePereiro/Chemostat_Heerden2013.jl in the Julia Pkg REPL to install the
    # package, then you must activate the package enviroment (see README)
    import Chemostat_Heerden2013
    const ChH = Chemostat_Heerden2013
    const Hd  = ChH.HeerdenData;
    const Bd  = ChH.BegData
    const iJR = ChH.iJR904
end

## -------------------------------------------------------------------
# BASE MODEL
# Load Mat file
function load_src_model()
    src_file = iJR.rawdir("iJR904.mat")
    mat_model = MAT.matread(src_file)["model"]
   return ChU.MetNet(mat_model; reshape = true);
end

get_comb(n::Int) = (digits(i; base = 2, pad = n) for i in ((2^n) - 1):-1:0)

## -------------------------------------------------------------------
# globals
inac_dat_cid = ("INACT_DAT")
inac_map = iJR.load_Hd_to_inactivate_map()
inac_rxns = vcat(values(inac_map)...)

## -------------------------------------------------------------------
# Collect data
let
    # UJL.delete_cache(inac_dat_cid) # Uncomment to recompute
    UJL.exist_cache(inac_dat_cid) && return

    inac_comb_pool = get_comb(length(inac_rxns))
    
    # init val
    # model0 = iJR.load_model("base_model")
    model0 = load_src_model()
    biom_ider = iJR.BIOMASS_IDER
    fbaout = ChLP.fba(model0, biom_ider)
    biom_val0 = ChU.av(model0, fbaout, biom_ider)
    inac_idxs = map(inac_rxns) do rxn
        ChU.rxnindex(model0, rxn)
    end

    model_pools = deepcopy(model0)
    prog = ProgressUnknown(;dt = 0.5)
    biom_vals = []
    for (ci, inac_comb) in enumerate(inac_comb_pool)

        to_inac = inac_rxns[isone.(inac_comb)]
        
        # inactivate
        for rxn in to_inac
            ChU.bounds!(model, rxn, 0.0, 0.0)
        end

        # Test FBA
        try
            fbaout = ChLP.fba(model, biom_ider)
            biom_val = ChU.av(model, fbaout, biom_ider)
            next!(prog; showvalues = [ ("biom_val", biom_val)])    
            push!(biom_vals, biom_val)
        catch err;
            next!(prog; showvalues = [("biom_val", "Fails")])     
            push!(biom_vals, 0.0)
        end

        # restart model
        model.lb[inac_idxs] .= model0.lb[inac_idxs]
        model.ub[inac_idxs] .= model0.ub[inac_idxs]

    end

    UJL.save_cache(inac_dat_cid, biom_vals)
end

## -------------------------------------------------------------------
# Plots
# growth vs n_inac
let
    biom_vals = UJL.load_cache(inac_dat_cid; verbose = false);
    inac_comb_pool = get_comb(length(inac_rxns))

    growth_hist = Dict()
    for (comb, biom_val) in zip(inac_comb_pool, biom_vals)
        ninac = sum(comb)
        get!(growth_hist, ninac, 0.0)
        growth_hist[ninac] = max(growth_hist[ninac], biom_val)
    end

    # plot(sum.(inac_comb_pool), biom_vals; label = "")
    bar(growth_hist; xlabel = "ninac", ylabel = "max growth", label = "")
    
end

## -------------------------------------------------------------------
# Plots
# growth vs inac rxn
let
    biom_vals = UJL.load_cache(inac_dat_cid; verbose = false);
    inac_comb_pool = get_comb(length(inac_rxns))

    growth_hist = Dict()
    for (comb, biom_val) in zip(inac_comb_pool, biom_vals)
        comb_inac_rxns = inac_rxns[isone.(comb)]
        for rxn in comb_inac_rxns
            get!(growth_hist, rxn, 0.0)
            growth_hist[rxn] += biom_val
        end
    end

    # plot(sum.(inac_comb_pool), biom_vals; label = "")
    bar(growth_hist; 
        xlabel = "inac rxn", ylabel = "accum. growth", 
        label = "", 
        xrotation = 35
    )
    
end