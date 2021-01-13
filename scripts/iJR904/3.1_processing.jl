import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

# Run add https://github.com/josePereiro/Chemostat_Heerden2013.jl in the Julia Pkg REPL to install the
# package, then you must activate the package enviroment (see README)
import Chemostat_Heerden2013
const Hd  = Chemostat_Heerden2013.HeerdenData
const BD  = Chemostat_Heerden2013.BegData
const iJR = Chemostat_Heerden2013.iJR904

# run add "https://github.com/josePereiro/Chemostat" in the 
# julia Pkg REPL for installing the package
import Chemostat
import Chemostat.LP.MathProgBase
const ChU = Chemostat.Utils
const ChSS = Chemostat.SteadyState
const ChLP = Chemostat.LP
const ChEP = Chemostat.MaxEntEP
const ChSU = Chemostat.SimulationUtils
const ChP = Chemostat.Plots

import SparseArrays
using Plots
using Statistics

## -------------------------------------------------------------------
# save results
DATA = ChU.load_data(iJR.MAXENT_FBA_EB_BOUNDLES_FILE)
SDATA = sort(collect(DATA); by = first);
@info "Loaded" length(SDATA)

## -------------------------------------------------------------------
const fileid = "3.1"

## -------------------------------------------------------------------
mysavefig(p, fname::String) = (savefig(p, joinpath(iJR.MODEL_FIGURES_DIR, fname)); p)
mysavefig(fname::String, p) = mysavefig(p, fname)
Base.minmax(a::Vector) = (minimum(a), maximum(a))

## -------------------------------------------------------------------
# Collect data
CONC_IDERS = String["GLC", "AcA", "FA"]
FLX_IDERS = [CONC_IDERS; "BiomassEcoli"]
color_pool = Dict(
    "GLC" => :red, "SA" => :yellow, "AcA" => :pink, 
    "FA" => :orange, "MA" => :blue, "BiomassEcoli" => :black
)
P = Dict()
D = ChU.DictTree()
objider = iJR.BIOMASS_IDER

let 
    for (exp, D0) in SDATA

        model = D0["model"]
        objidx = ChU.rxnindex(model, objider)
        exp_beta = D0["exp_beta"]
        epout = D0["epouts"][exp_beta]
        exp_xi = Hd.val("xi", exp)
        fbaout = D0["fbaout"]

        # Biomass
        fba_biom = ChU.av(model, fbaout, objidx)
        ep_biom = ChU.av(model, epout, objidx)
        ep_std = sqrt(ChU.va(model, epout, objidx))
        Hd_biom = Hd.val("D", exp)
        fva_range = getindex.(D0["fva"], objidx)
        
        # store
        D[:fba  , :flx, objider, exp] = fba_biom
        D[:ep   , :flx, objider, exp] = ep_biom
        D[:eperr, :flx, objider, exp] = ep_std
        D[:Hd   , :flx, objider, exp] = Hd_biom
        D[:fva  , :flx, objider, exp] = fva_range
        
        D[:color, objider, exp] = color_pool[objider]
        
        # mets
        for Hd_met in Hd.MSD_METS

            try
                model_met = iJR.Hd_mets_map[Hd_met]
                model_exch = iJR.exch_met_map[model_met]
                model_exchi = ChU.rxnindex(model, model_exch)

                # fuxes
                fba_av = ChU.av(model, fbaout, model_exchi)
                ep_av = ChU.av(model, epout, model_exchi)
                ep_std = sqrt(ChU.va(model, epout, model_exchi))
                Hd_flx = Hd.val("u$Hd_met", exp)
                fva_range = getindex.(D0["fva"], model_exchi)
                
                # conc (s = c + u*xi)
                c = Hd.val("c$Hd_met", exp, 0.0)
                fba_conc = max(c + fba_av * exp_xi, 0.0)
                ep_conc = max(c + ep_av * exp_xi, 0.0)
                Hd_conc = Hd.val("s$Hd_met", exp)
                
                for (datk, datt, dat) in [(:fba, :flx, fba_av), (:fba, :conc, fba_conc), 
                                        (:ep, :flx, ep_av), (:ep, :conc, ep_conc), 
                                        (:eperr, :flx, ep_std), (:eperr, :conc, ep_std * exp_xi), 
                                        (:Hd, :flx, Hd_flx), (:Hd, :conc, Hd_conc), (:fva, :flx, fva_range)]
                    D[datk, datt, Hd_met, exp] = dat
                end
                D[:color, Hd_met, exp] = color_pool[Hd_met]
            
            catch err
                @warn string(Hd_met, " fails") err
            end
        end

    end

    D[:mflx], D[:Mflx] = minmax(D[[:ep, :fba, :Hd], :flx, FLX_IDERS, Hd.EXPS])
    D[:mconc], D[:Mconc] = minmax(D[[:ep, :fba, :Hd], :conc, CONC_IDERS, Hd.EXPS])

end

## -------------------------------------------------------------------
# beta_exp corr
let
    p = plot(title = iJR.PROJ_IDER, xlabel = "exp biom", ylabel = "model biom")
    ep_vals = D[:ep, :flx, iJR.BIOMASS_IDER, Hd.EXPS]
    Hd_vals = D[:Hd, :flx, iJR.BIOMASS_IDER, Hd.EXPS]
    color = D[:color, iJR.BIOMASS_IDER, Hd.EXPS]
    scatter!(p, ep_vals, Hd_vals; label = "", color)
    P["obj_val_ep_corr"] = deepcopy(p)
    mysavefig(p, "$(fileid)_obj_val_ep_corr.png")
end

## -------------------------------------------------------------------
# flux vs beta
let
    p = plot(title = iJR.PROJ_IDER, xlabel = "beta", ylabel = "biom")
    for (exp, D0) in SDATA
        model = D0["model"]
        objidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
        epouts = D0["epouts"]
        exp_beta = D0["exp_beta"]
        exp_xi = D0["exp_xi"]
        scatter!(p, [exp_beta], [Hd.val("D", exp)], ms = 12, color = :white, label = "")

        betas = collect(keys(epouts)) |> sort
        bioms = [ChU.av(model, epouts[beta], objidx) for beta in betas]
        scatter!(p, betas, bioms, label = "", color = :black, alpha = 0.2)

    end
    P["obj_val_vs_beta"] = deepcopy(p)
    mysavefig(p, "$(fileid)_obj_val_vs_beta.png")
end


## -------------------------------------------------------------------
# total correlations
let
    for (dat_prefix, iders, zoom_lim) in [(:flx, FLX_IDERS, [-2.5, 2.5]), 
                                            (:conc, CONC_IDERS, [0.0, 100.0])]

        ep_vals = D[:ep, dat_prefix, iders, Hd.EXPS]
        ep_errs = D[:eperr, dat_prefix, iders, Hd.EXPS]
        fba_vals = D[:fba, dat_prefix, iders, Hd.EXPS]
        Hd_vals = D[:Hd, dat_prefix, iders, Hd.EXPS]
        color = D[:color, iders, Hd.EXPS]
        m, M = minmax([fba_vals; ep_vals; Hd_vals])

        scatter_params = (;label = "", color, ms = 8, alpha = 0.8)
        # ep corr
        p1 = plot(title = "$(iJR.PROJ_IDER) (EP)", 
            ylabel = "exp $(dat_prefix)",
            xlabel = "model $(dat_prefix)", 
        )
        scatter!(p1, ep_vals, Hd_vals; xerr = ep_errs, scatter_params...)
        plot!(p1, [m,M], [m,M]; ls = :dash, color = :black, label = "")

        # fba corr
        p2 = plot(title = "$(iJR.PROJ_IDER) (FBA)", 
            ylabel = "exp $(dat_prefix)",
            xlabel = "model $(dat_prefix)", 
        )
        scatter!(p2, fba_vals, Hd_vals; scatter_params...)
        plot!(p2, [m,M], [m,M]; ls = :dash, color = :black, label = "")

        p = plot(p1, p2; layout = 2)
        pname = string(dat_prefix, "_tot_corr")
        P[pname] = deepcopy(p)
        mysavefig(p, "$(fileid)_$(pname).png")

        # zoom
        p = plot(p1, p2; layout = 2, xlim = zoom_lim, ylim = zoom_lim)
        pname = string(dat_prefix, "_tot_corr_zoom")
        P[pname] = deepcopy(p)
        mysavefig(p, "$(fileid)_$(pname).png")
    end

end

## -------------------------------------------------------------------
# fva bounds
let
    ps = []
    for ider = FLX_IDERS
        p = plot(title = ider, xlabel = "replica", ylabel = "flx")
        Hd_vals = D[:Hd, :flx, ider, Hd.EXPS]
        ep_vals = D[:ep, :flx, ider, Hd.EXPS]
        fba_vals = D[:fba, :flx, ider, Hd.EXPS]
        fva_ranges = D[:fva, :flx, ider, Hd.EXPS]
        xticks =  (Hd.EXPS, string.(Hd.EXPS))
        plot!(p, Hd.EXPS, first.(fva_ranges); 
            label = "lb", color = :blue, alpha = 0.8, ls = :dot, lw = 3, xticks)
        plot!(p, Hd.EXPS, Hd_vals; 
            label = "exp", color = :black, alpha = 0.5, lw = 3, xticks)
        plot!(p, Hd.EXPS, ep_vals; 
            label = "ep", color = :blue, alpha = 0.5, lw = 5, ls = :dash, xticks)
        plot!(p, Hd.EXPS, fba_vals; 
            label = "fba", color = :red, alpha = 0.5, lw = 5, ls = :dash, xticks)
        plot!(p, Hd.EXPS, last.(fva_ranges);  
            label = "ub", color = :red, alpha = 0.8, ls = :dot, lw = 3, xticks)
        
        pname = string(ider, "_bound_study")
        P[pname] = deepcopy(p)
        mysavefig(p, "$(fileid)_$(pname).png")
        push!(ps, p)
    end
    plot(ps..., layout = length(ps))
end

## -------------------------------------------------------------------
# marginal distributions
let
    
    for Hd_met in CONC_IDERS
        p0 = plot(xaxis = nothing, yaxis = nothing, titlefont = 10, plot_title = Hd_met)
        ps = []
        for (exp, D0) in SDATA
            p = deepcopy(p0)
            model = D0["model"]
            objidx = ChU.rxnindex(model, objider)
            exp_beta = D0["exp_beta"]
            epout = D0["epouts"][exp_beta]
            exp_xi = Hd.val("xi", exp)
            fbaout = D0["fbaout"]

            model_met = iJR.Hd_mets_map[Hd_met]
            model_exch = iJR.exch_met_map[model_met]
                    
            ChP.plot_marginal!(p, model, [epout, fbaout], model_exch; 
                title = exp, legend = false)
            vline!(p, [Hd.val("u$(Hd_met)", exp)]; label = "", lw = 6, color = :black)

            m, M = minmax(D[[:ep, :fba, :Hd], :flx, Hd_met, exp])
            margin = 3 * D[:eperr, :flx, Hd_met, exp]
            plot!(p; xlim = [m - margin, M + margin])
            push!(ps, p)
        end

        for k in [:xi, :D, :cGLC]
            p = deepcopy(p0)
            p = plot!(p, Hd.EXPS, Hd.val(k); title = k, label = "", lw = 5, ylim = [0.0, Inf])
            push!(ps, p)
        end
        
        p = plot(ps...; layout = length(ps))

        pname = string(Hd_met, "_marginals")
        P[pname] = deepcopy(p)
        mysavefig(p, "$(fileid)_$(pname).png")

    end

end 