import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid

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
    import JuMP.MathOptInterface
    using Serialization
    import UtilsJL
    const UJL = UtilsJL
    
    import FileIO
    using Plots
    import GR
    GR.inline("png")
end

## -----------------------------------------------------------------------------------------------
LPDAT = ChU.load_data(iJR.LP_DAT_FILE)
const FBA_BOUNDED = :FBA_BOUNDED
const FBA_OPEN = :FBA_OPEN
const YIELD = :YIELD

## -----------------------------------------------------------------------------------------------
fileid = "3.1"
function mysavefig(p, pname; params...) 
    fname = UJL.mysavefig(p, string(fileid, "_", pname), iJR.MODEL_FIGURES_DIR; params...)
    @info "Plotting" fname
end
myminmax(a::Vector) = (minimum(a), maximum(a))  
CONC_IDERS = ["GLC", "SA", "AcA", "FA"]
FLX_IDERS = ["GLC", "SA", "AcA", "FA"]

EXPS = Hd.EXPS 

exp_colors = let
    colors = Plots.distinguishable_colors(length(EXPS))
    Dict(exp => color for (exp, color) in zip(EXPS, colors))
end

ider_colors = Dict(
    "GLC" => :red, "SA" => :yellow,
    "AcA" => :orange, "FA" => :blue,
    "D" => :black,
)

method_colors = Dict(
    FBA_OPEN => :red,
    FBA_BOUNDED => :orange,
    YIELD => :blue,
)

## -----------------------------------------------------------------------------------------------
# yield correlation
let
    p = plot(;title = "Yield correlation", xlabel = "exp", ylabel = "model")
    m, M = Inf, -Inf
    for exp in EXPS
        !haskey(LPDAT, YIELD, :model, exp) && continue
        model = LPDAT[YIELD, :model, exp]
        yout = LPDAT[YIELD, :yout, exp]
        model_yield = LPDAT[YIELD, :yield, exp]
        
        objidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
        exp_growth = Hd.val("D", exp)
        model_growth = ChU.av(model, yout, objidx)

        diff = abs(model_growth - exp_growth)/exp_growth
        diff > 0.05 && continue # unfeasible

        exp_yield = abs(Hd.val("D", exp) / Hd.uval("GLC", exp))
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
    for id in [:D, :cGLC, :xi, :uGLC, :sGLC]
        p = plot(;title = "yield vs $(id)", xlabel = "exp $id", ylabel = "yield")
        for exp in EXPS
            !haskey(LPDAT, YIELD, :model, exp) && continue
            model = LPDAT[YIELD, :model, exp]
            yout = LPDAT[YIELD, :yout, exp]
            model_yield = LPDAT[YIELD, :yield, exp]
            # status, yflxs, model_yield, d, model = DAT[D]
            exglc_idx = ChU.rxnindex(model, "EX_glc_LPAREN_e_RPAREN__REV")
            biomass_idx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
            
            Hd_val = Hd.val(id, exp)
            exp_yield = abs(Hd.val("D", exp) / Hd.uval("GLC", exp))
            scatter!(p, [Hd_val], [model_yield]; ms = 8,
                    color = :blue, alpha = 0.6, label = ""
            )
            scatter!(p, [Hd_val], [exp_yield]; ms = 8,
                    color = :red, alpha = 0.6, label = ""
            )
        end
        push!(ps, p)
    end
    pname = "yield_vs_stuff"
    mysavefig(ps, pname)
end

## -----------------------------------------------------------------------------------------------
# correlations
DAT = UJL.DictTree()
DAT[:FLX_IDERS] = FLX_IDERS
DAT[:CONC_IDERS] = CONC_IDERS
DAT[:EXPS] = EXPS

## -----------------------------------------------------------------------------------------------
# ["GLC", "AcA", "FA"]
FLX_IDERS_MAP = Dict(
    "GLC" => "EX_glc_LPAREN_e_RPAREN__REV",
    "AcA" => "EX_ac_LPAREN_e_RPAREN_",
    "FA" => "EX_for_LPAREN_e_RPAREN_",
    "SA" => "EX_succ_LPAREN_e_RPAREN_",
)

# Flx correlations
let

    yield_p = plot(title = "yield tot corrs"; xlabel = "exp flx", ylabel = "model flx")
    open_fba_p = plot(title = "open fba tot corrs"; xlabel = "exp flx", ylabel = "model flx")
    bounded_fba_p = plot(title = "bounded fba tot corrs"; xlabel = "exp flx", ylabel = "model flx")

    margin, m, M = -Inf, Inf, -Inf
    for (Hd_ider, model_ider) in FLX_IDERS_MAP
        for exp in EXPS

                color = ider_colors[Hd_ider]
                Hd_flx = abs(Hd.uval(Hd_ider, exp)) # every body is possitive here
                
                # yield
                !haskey(LPDAT, YIELD, :model, exp) && continue
                model = LPDAT[YIELD, :model, exp]
                yout = LPDAT[YIELD, :yout, exp]

                ymax_flx = ChU.av(model, yout, model_ider)
                DAT[YIELD, :Hd, :flx, Hd_ider, exp] = Hd_flx
                DAT[YIELD, :lp, :flx, Hd_ider, exp] = ymax_flx

                scatter!(yield_p, [Hd_flx], [ymax_flx]; ms = 8,
                    color, alpha = 0.6, label = ""
                )
                
                # bounded fba
                for (fba_type, p) in [(FBA_BOUNDED, bounded_fba_p) , 
                                    (FBA_OPEN, open_fba_p)]

                    model = LPDAT[fba_type, :model, exp]
                    fbaout = LPDAT[fba_type, :fbaout, exp]
                    
                    fba_flx = ChU.av(model, fbaout, model_ider)
                    DAT[fba_type, :Hd, :flx, Hd_ider, exp] = Hd_flx
                    DAT[fba_type, :lp, :flx, Hd_ider, exp] = fba_flx

                    scatter!(p, [Hd_flx], [fba_flx]; ms = 8,
                        color, alpha = 0.6, label = ""
                    )
                    m = minimum([m, Hd_flx, ymax_flx, fba_flx])
                    M = maximum([M, Hd_flx, ymax_flx, fba_flx])
                end

        end
    end
    margin = abs(M - m) * 0.1
    ps = [yield_p, bounded_fba_p, open_fba_p]
    for p in ps
        plot!(p, [m - margin, M + margin], [m - margin, M + margin]; 
            ls = :dash, color = :black, label = "")
    end
    
    pname = "flx_tot_corr"
    layout = (1, 3)
    mysavefig(ps, pname; layout)

end

## -----------------------------------------------------------------------------------------------
# Dev
let
    dsource = :Hd
    msource = :lp
    method = FBA_BOUNDED
    D = DAT

    exp_vals, model_vals = [], []
    for exp in D[:EXPS], ider in D[:FLX_IDERS]
        if !haskey(D, method, msource, :flx, ider, exp) 
            @warn("Not found", ider, exp)
            continue
        end

        exp_val = D[method, dsource, :flx, ider, exp]
        model_val = D[method, msource, :flx, ider, exp]
        push!(exp_vals, exp_val)
        push!(model_vals, model_val)
    end

    diffsign = sign.(exp_vals) .* sign.(model_vals)
    diffsign = ifelse.(diffsign .== 0, 1.0, diffsign)
    exp_vals = abs.(exp_vals) .* diffsign
    model_vals = abs.(model_vals) .* diffsign
    
    p = scatter(exp_vals, model_vals; label = "false")
    plot!(p, exp_vals, exp_vals; label = "", ls = :dash)
    pname = "test.png" 
    mysavefig(p, pname)

end

## -------------------------------------------------------------------
# Conc correlations
let
    yield_p = plot(title = "yield tot corrs"; xlabel = "exp conc", ylabel = "model conc")
    open_fba_p = plot(title = "open fba tot corrs"; xlabel = "exp conc", ylabel = "model conc")
    bounded_fba_p = plot(title = "bounded fba tot corrs"; xlabel = "exp conc", ylabel = "model conc")
    
    margin, m, M = -Inf, Inf, -Inf
    for (Hd_ider, model_ider) in FLX_IDERS_MAP
        for exp in EXPS

                color = ider_colors[Hd_ider]
                Hd_sval = Hd.sval(Hd_ider, exp)
                Hd_cval = Hd.cval(Hd_ider, exp, 0.0)
                exp_xi = Hd.val(:xi, exp)

                # yield
                !haskey(LPDAT, YIELD, :model, exp) && continue
                model = LPDAT[YIELD, :model, exp]
                yout = LPDAT[YIELD, :yout, exp]

                # conc (s = c + u*xi)
                ymax_flx = ChU.av(model, yout, model_ider)
                ymax_sval = Hd_cval == 0 ? ymax_flx * exp_xi :
                    max(Hd_cval - ymax_flx * exp_xi, 0.0)
                DAT[YIELD, :Hd, :conc, Hd_ider, exp] = Hd_sval
                DAT[YIELD, :lp, :conc, Hd_ider, exp] = ymax_flx    

                scatter!(yield_p, [Hd_sval], [ymax_sval]; ms = 8,
                    color, alpha = 0.6, label = ""
                )

                # fba
                for (fba_type, p) in [(FBA_BOUNDED, bounded_fba_p) , 
                                    (FBA_OPEN, open_fba_p)]

                    model = LPDAT[fba_type, :model, exp]
                    fbaout = LPDAT[fba_type, :fbaout, exp]
                    
                    # conc (s = c + u*xi)
                    fba_flx = ChU.av(model, fbaout, model_ider)
                    fba_sval = Hd_cval == 0 ? fba_flx * exp_xi :
                        max(Hd_cval - fba_flx * exp_xi, 0.0)
                    DAT[fba_type, :Hd, :conc, Hd_ider, exp] = Hd_sval
                    DAT[fba_type, :lp, :conc, Hd_ider, exp] = fba_sval

                    scatter!(p, [Hd_sval], [fba_sval]; ms = 8,
                        color, alpha = 0.6, label = ""
                    )
                    m = minimum([m, Hd_sval, ymax_sval, fba_sval])
                    M = maximum([M, Hd_sval, ymax_sval, fba_sval])
                end

        end
    end

    margin = abs(M - m) * 0.1
    ps = [yield_p, bounded_fba_p, open_fba_p]
    for p in ps
        plot!(p, [m - margin, M + margin], [m - margin, M + margin]; 
            ls = :dash, color = :black, label = "")
    end
    
    pname = "conc_tot_corr"
    layout = (1, 3)
    mysavefig(ps, pname; layout)
end

## -------------------------------------------------------------------
# Inter project comunication
let
    CORR_DAT = isfile(iJR.CORR_DAT_FILE) ? ChU.load_data(iJR.CORR_DAT_FILE) : Dict()
    CORR_DAT[:LP] = DAT
    ChU.save_data(iJR.CORR_DAT_FILE, CORR_DAT)
end

## -------------------------------------------------------------------
# join flx correlations
let

    figdir = iJR.MODEL_FIGURES_DIR
    for (lp_p, ep_p, join_name) in [
        ("3.1_flx_tot_corr.png", "2.1_flx_tot_corr.png", "flx_join_corr.png"),
        ("3.1_conc_tot_corr.png", "2.1_conc_tot_corr.png", "conc_join_corr.png"),
    ] 
        lp_img = FileIO.load(joinpath(figdir, lp_p))
        ep_img = FileIO.load(joinpath(figdir, ep_p))
        join_p = UJL.make_grid([lp_img, ep_img])
        fname = joinpath(figdir, join_name)
        FileIO.save(fname, join_p)
        @info "Plotting" fname
    end
end

## -------------------------------------------------------------------
# # # leyends
# # # TODO fix this...
# # let
# #     for (title, colors) in [
# #             ("exp", exp_colors), 
# #             ("iders", ider_colors),
# #             ("method", method_colors)
# #         ]
# #     p = plot(; framestyle = :none)
# #         scolors = sort(collect(colors); by = (p) -> string(first(p)))
# #         for (id, color) in scolors
# #             scatter!(p, [0], [0];
# #                 thickness_scaling = 1,
# #                 color, ms = 8, label = string(id),
# #                 legendfontsize=10, 
# #                 # size = [300, 900],
# #                 # legend = :left
# #             )
# #         end
# #         mysavefig(p, "$(title)_color_legend")
# #     end
# # end