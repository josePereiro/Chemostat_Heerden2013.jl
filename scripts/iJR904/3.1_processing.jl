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

const FBA_MIN_COST = :FBA_MIN_COST
const FBA_OPEN = :FBA_OPEN
const YIELD = :YIELD
const FBA_MAX_VG_YIELD = :FBA_MAX_VG_YIELD

## -----------------------------------------------------------------------------------------------
fileid = "3.1"
function mysavefig(p, pname; params...) 
    fname = UJL.mysavefig(p, string(fileid, "_", pname), iJR.MODEL_FIGURES_DIR; params...)
    @info "Plotting" fname
end
myminmax(a::Vector) = (minimum(a), maximum(a))  
CONC_IDERS = ["GLC", "SUCC", "AC", "FORM"]
FLX_IDERS = ["GLC", "SUCC", "AC", "FORM"]

EXPS = Hd.EXPS 

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
    FBA_OPEN => :red,
    FBA_MIN_COST => :orange,
    YIELD => :blue,
)

ALL_METHODS = [
    FBA_MAX_VG_YIELD,
    FBA_MIN_COST, 
    FBA_OPEN
]

Hd_mets_map = iJR.load_Hd_mets_map()
Hd_rxns_map = iJR.load_Hd_rxns_map()

## -----------------------------------------------------------------------------------------------
# correlations
DAT = UJL.DictTree()
DAT[:FLX_IDERS] = FLX_IDERS
DAT[:CONC_IDERS] = CONC_IDERS
DAT[:EXPS] = EXPS

## -----------------------------------------------------------------------------------------------
# Flx correlations
let
    ps = Plots.Plot[]
    for method in ALL_METHODS
        p = plot(;title = string(method), xlabel = "exp flx", ylabel = "model flx")
        for Hd_ider in FLX_IDERS
            model_ider = Hd_rxns_map[Hd_ider]
            for exp in EXPS
                color = ider_colors[Hd_ider]
                # every body is possitive here
                Hd_flx = abs(Hd.uval(Hd_ider, exp))

                model = LPDAT[method, :model, exp]
                fbaout = LPDAT[method, :fbaout, exp]
                
                fba_flx = abs(ChU.av(model, fbaout, model_ider))
                DAT[method, :Hd, :flx, Hd_ider, exp] = Hd_flx
                DAT[method, :lp, :flx, Hd_ider, exp] = fba_flx

                scatter!(p, [Hd_flx], [fba_flx]; 
                    ms = 8, color, label = ""
                )
            end
        end
        xs = DAT[method, [:Hd, :lp], :flx, FLX_IDERS, EXPS] |> sort
        plot!(p, xs, xs; label = "", ls = :dash, 
            alpha = 0.9, color = :black, lw = 3
        )
        push!(ps, p)
    end
    
    pname = "flx_tot_corr"
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