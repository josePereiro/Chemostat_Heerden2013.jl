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
    using Plots
end

## -----------------------------------------------------------------------------------------------
INDEX = ChU.load_data(iJR.MAXENT_VARIANTS_INDEX_FILE; verbose = false);
const HOMO = :HOMO
const BOUNDED = :BOUNDED
const EXPECTED = :EXPECTED

fileid = "5"
function mysavefig(p, pname; params...) 
    fname = UJL.mysavefig(p, string(fileid, "_", pname), iJR.MODEL_FIGURES_DIR; params...)
    @info "Plotting" fname
end
myminmax(a::Vector) = (minimum(a), maximum(a))
CONC_IDERS = String["GLC", "AcA", "FA"]
FLX_IDERS = [CONC_IDERS; "D"]

FLX_IDERS_MAP = Dict(
    "uGLC" => "EX_glc_LPAREN_e_RPAREN__REV",
    "uAcA" => "EX_ac_LPAREN_e_RPAREN_",
    "uFA" => "EX_for_LPAREN_e_RPAREN_",
    "D" => "BiomassEcoli"
)

exp_colors = let
    cGLCs = Hd.val("cGLC")
    colors = Plots.distinguishable_colors(length(cGLCs))
    Dict(exp => color for (exp, color) in zip(eachindex(cGLCs), colors))
end

ider_colors = let
    colors = Plots.distinguishable_colors(length(FLX_IDERS))
    Dict(met => color for (met, color) in zip(FLX_IDERS, colors))
end

method_colors = Dict(
    HOMO => :red,
    BOUNDED => :orange,
    EXPECTED => :blue,
)

## -----------------------------------------------------------------------------------------------
function yLP(S, b, lb, ub, c, d; 
        ϵ = 1e-5, 
        sense = JuMP.MOI.MAX_SENSE,
        solver = GLPK.Optimizer
    )

    lp_model = JuMP.Model(solver)
    M, N = size(S)

    # Variables
    y = JuMP.@variable(lp_model, y[1:N])
    r = JuMP.@variable(lp_model, r)

    # Constraints
    JuMP.@constraint(lp_model, S * y - r .* b .== 0.0)
    JuMP.@constraint(lp_model, d' * y == 1.0)
    JuMP.@constraint(lp_model, y - r .* lb .>= 0.0)
    JuMP.@constraint(lp_model, y - r .* ub .<= 0.0)
    JuMP.@constraint(lp_model, r >= ϵ)

    # objective
    JuMP.@objective(lp_model, sense, c' * y)

    # optimize
    JuMP.optimize!(lp_model)

    yval = JuMP.value.(y)
    rval = JuMP.value.(r)
    v = yval ./ rval # sol
    y = (c' * v) / (d' * v) # yield

    status = JuMP.termination_status(lp_model)
    status, v, y
end
yLP(net, d; kwargs...) = yLP(net.S, net.b, net.lb, net.ub, net.c, d; kwargs...)

## -----------------------------------------------------------------------------------------------
# BASE_MODELS = ChU.load_data(iJR.BASE_MODELS_FILE; verbose = false);
# model_dict = BASE_MODELS["base_model"]
# model = ChU.MetNet(;model_dict...) |> ChU.uncompressed_model

# for method in [HOMO, EXPECTED, BOUNDED]
DAT = UJL.DictTree()
let
    objider = iJR.BIOMASS_IDER
    method = HOMO
    for exp in Hd.EXPS
        try
            # load
            datfile = INDEX[method, :DFILE, exp]
            dat = deserialize(datfile)
            
            # prepare model
            model = dat[:model] |> ChU.uncompressed_model
            M, N = size(model)
            model = ChU.MetNet(model; c = zeros(N))
            ChU.check_dims(model)
            ChU.invert_bkwds!(model)
            exglc_idx = ChU.rxnindex(model, "EX_glc_LPAREN_e_RPAREN__REV")
            exac_idx = ChU.rxnindex(model, "EX_ac_LPAREN_e_RPAREN_")
            biomass_idx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
            exp_growth = Hd.val("D", exp)
            dgrowth = 0.0
            # ChU.lb!(model, iJR.BIOMASS_IDER, exp_growth * (1 - dgrowth))
            ChU.ub!(model, iJR.BIOMASS_IDER, exp_growth * (1 + dgrowth))

            # fba
            fbaout = ChLP.fba(model, iJR.BIOMASS_IDER)
            fba_growth = ChU.av(model, fbaout, iJR.BIOMASS_IDER)

            # yield max
            model.c[biomass_idx] = 1.0
            d = zeros(N); 
            d[exglc_idx] = 1.0
            # d[exac_idx] = 1.0
            status, yflxs, yield = yLP(model, d)
            status != JuMP.MOI.OPTIMAL && @warn status
            ymax_growth = yflxs[biomass_idx]
            DAT[exp] = (;status, yflxs, yield, d, model)

            @info("Yield Maximization", 
                exp, status, yield,
                fba_growth, ymax_growth, exp_growth
            ); println()

        catch err
            @warn("Error", err, exp); println()
        end
    end
end


## -----------------------------------------------------------------------------------------------
# yield correlation
let
    p = plot(;title = "Yield correlation", xlabel = "exp", ylabel = "model")
    m, M = Inf, -Inf
    for exp in Hd.EXPS
        !haskey(DAT, exp) && continue
        status, yflxs, model_yield, d, model = DAT[exp]
        exglc_idx = ChU.rxnindex(model, FLX_IDERS_MAP["uGLC"])
        biomass_idx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
        
        exp_yield = abs(Hd.val("D", exp) / Hd.val("uGLC", exp))
        scatter!(p, [exp_yield], [model_yield]; ms = 8,
                color = :blue, alpha = 0.6, label = ""
        )
        m = minimum([m, exp_yield, model_yield])
        M = maximum([M, exp_yield, model_yield])
    end
    plot!(p, [m,M], [m,M]; ls = :dash, color = :black, label = "")
    pname = "yield_corr"
    mysavefig(p, pname)
end

## -----------------------------------------------------------------------------------------------
# yield vs stuff
let
    ps = Plots.Plot[]
    for id in ["D", "cGLC", "xi", "Ysx", "uGLC"]
        p = plot(;title = "Yield vs $(id)", xlabel = id, ylabel = "yield")
        for exp in Hd.EXPS
            !haskey(DAT, exp) && continue
            status, yflxs, model_yield, d, model = DAT[exp]
            exglc_idx = ChU.rxnindex(model, "EX_glc_LPAREN_e_RPAREN__REV")
            biomass_idx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
            
            Hd_val = Hd.val(id, exp)
            scatter!(p, [Hd_val], [model_yield]; ms = 8,
                    color = :blue, alpha = 0.6, label = ""
            )
        end
        push!(ps, p)
    end
    pname = "yield_vs_stuff"
    mysavefig(ps, pname)
end
## -----------------------------------------------------------------------------------------------
# correlations
let

    p = plot(title = "tot corrs"; xlabel = "exp", ylabel = "model")

    for (Hd_ider, model_ider) in FLX_IDERS_MAP
        margin, m, M = -Inf, Inf, -Inf
        for exp in Hd.EXPS
            !haskey(DAT, exp) && continue
            status, yflxs, yield, d, model = DAT[exp]
            model_idx = ChU.rxnindex(model, model_ider)
            
            sense = Hd_ider == "D" ? 1 : -1
            Hd_flx = sense * Hd.val(Hd_ider, exp)
            ymax_flx = yflxs[model_idx]
            scatter!(p, [Hd_flx], [ymax_flx]; 
                color = :blue, alpha = 0.6, label = ""
            )
        end
    end

    pname = "tot_corr"
    mysavefig(p, pname)
end