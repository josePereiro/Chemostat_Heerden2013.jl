import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

# Run add https://github.com/josePereiro/Chemostat_Heerden2013.jl in the Julia Pkg REPL to install the
# package, then you must activate the package enviroment (see README)
import Chemostat_Heerden2013
const Hd  = Chemostat_Heerden2013.HeerdenData;
const BD  = Chemostat_Heerden2013.BegData;
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

import SparseArrays
using Plots

## -------------------------------------------------------------------
# save results
DATA = ChU.load_data(iJR.MAXENT_FBA_EB_BOUNDLES_FILE)
SDATA = sort(collect(DATA); by = first);
@info "Loaded" length(SDATA)

## -------------------------------------------------------------------
const fileid = replace(basename(@__FILE__), ".jl" => "")

## -------------------------------------------------------------------
mysavefig(p, fname::String) = (savefig(p, joinpath(iJR.MODEL_FIGURES_DIR, fname)); p)
mysavefig(fname::String, p) = mysavefig(p, fname)

## -------------------------------------------------------------------
# Collect data
D = Dict()
objider = iJR.BIOMASS_IDER

function getpush!(k1::Symbol, kvs::Tuple...) 
    main_dict = get!(D, k1, Dict())
    for (k2, v) in kvs
        a = get!(main_dict, k2, [])
        push!(a, v)
    end
end

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
        
        # store
        metflxkey = (:flx, objider)
        expflxkey = (:flx, exp)
        getpush!(:fba, (expflxkey, fba_biom), (metflxkey, fba_biom))
        getpush!(:ep, (expflxkey, ep_biom), (metflxkey, ep_biom))
        getpush!(:eperr, (expflxkey, ep_std), (metflxkey, ep_std))
        getpush!(:Hd, (expflxkey, Hd_biom), (metflxkey, Hd_biom))
        
        # mets
        for Hd_met in Hd.msd_mets
            try
                model_met = iJR.Hd_mets_map[Hd_met]
                model_exch = iJR.exch_met_map[model_met]

                # fuxes
                fba_av = ChU.av(model, fbaout, model_exch)
                ep_av = ChU.av(model, epout, model_exch)
                ep_std = sqrt(ChU.va(model, epout, model_exch))
                Hd_flx = Hd.val("u$Hd_met", exp)
                
                # conc (s = c + u*xi)
                c = Hd.val("c$Hd_met", exp, 0.0)
                fba_conc = max(c + fba_av * exp_xi, 0.0)
                ep_conc = max(c + ep_av * exp_xi, 0.0)
                Hd_conc = Hd.val("s$Hd_met", exp)
                
                # store
                expflxkey = (:flx, exp)
                metflxkey = (:flx, Hd_met)
                allflxkey = (:flx, :all)
                metconckey = (:conc, Hd_met)
                expconckey = (:conc, exp)
                allconckey = (:conc, :all)

                for (k, flx, conc) in [
                                (:ep, ep_av, ep_conc), 
                                (:fba, fba_av, fba_conc),
                                (:eperr, ep_std, ep_std * exp_xi),
                                (:Hd, Hd_flx, Hd_conc)
                            ]
                        getpush!(k, 
                            (metflxkey, flx), (expflxkey, flx), (allflxkey, flx), 
                            (metconckey, conc), (expconckey, conc), (allconckey, conc), 
                        )
                end

                
                D[:mflx] = minimum([get!(D, :mflx, Inf), fba_av, ep_av, fba_av])
                D[:Mflx] = maximum([get!(D, :Mflx, -Inf), fba_av, ep_av, fba_av])
                D[:mconc] = minimum([get!(D, :mconc, Inf), fba_conc, ep_conc, Hd_conc])
                D[:Mconc] = maximum([get!(D, :Mconc, -Inf), fba_conc, ep_conc, Hd_conc])
            
            catch err
                @warn string(Hd_met, " fails") err
            end
        end

    end
end

## -------------------------------------------------------------------
# beta_exp corr
let
    p = plot(title = iJR.PROJ_IDER, xlabel = "exp biom", ylabel = "model biom")
    for (exp, D0) in SDATA

        model = D0["model"]
        objidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
        exp_beta = D0["exp_beta"]
        epout = D0["epouts"][exp_beta]

        ep_objval = ChU.av(model, epout, objidx)
        Hd_objval = Hd.val("D", exp)
        scatter!(p, [ep_objval], [Hd_objval], label = "")
    end
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
    mysavefig(p, "$(fileid)_obj_val_vs_beta.png")
end


## -------------------------------------------------------------------
# total correlations
let

    # plots
    function plot_res(;xlabel, ylabel, dat_prefix, norm_lims, lims, fname)

        ep_p = plot(;title = "$(iJR.PROJ_IDER) (EP)", xlabel, ylabel)
        ep_p_norm = plot(;title = "$(iJR.PROJ_IDER) (EP normalized)", xlabel, ylabel)
        fba_p = plot(;title = "$(iJR.PROJ_IDER) (FBA)", xlabel, ylabel)
        fba_p_norm = plot(;title = "$(iJR.PROJ_IDER) (FBA normalized)", xlabel, ylabel)

        for ider in [Hd.msd_mets; objider]
            # ider in ["SA", "MA"] && continue

            # flx
            flxkey = (dat_prefix, ider)
            if haskey(D[:Hd], flxkey)
            
                # fba
                scatter!(fba_p, D[:fba][flxkey], D[:Hd][flxkey]; label = "")
                norm = maximum([abs.(D[:fba][flxkey]); abs.(D[:Hd][flxkey])])
                scatter!(fba_p_norm, D[:fba][flxkey] ./ norm, D[:Hd][flxkey] ./ norm; label = "")

                # ep
                scatter!(ep_p, D[:ep][flxkey], D[:Hd][flxkey];
                    xerr = D[:eperr][flxkey], 
                    label = "")
                norm = maximum([abs.(D[:ep][flxkey]); abs.(D[:Hd][flxkey])])
                scatter!(ep_p_norm, D[:ep][flxkey] ./ norm, D[:Hd][flxkey] ./ norm; 
                    xerr = D[:eperr][flxkey] ./ norm, 
                    label = "")
            end
        end
        plot!(fba_p_norm, norm_lims, norm_lims; ls = :dash, color = :black, label = "")
        plot!(ep_p_norm, norm_lims, norm_lims; ls = :dash, color = :black, label = "")
        
        plot!(fba_p, lims, lims; ls = :dash, color = :black, label = "")
        plot!(ep_p, lims, lims; ls = :dash, color = :black, label = "")
        
        p = plot([ep_p_norm, fba_p_norm, ep_p, fba_p]...; layout = 4, size = [800, 800])
        mysavefig(p, fname)
    end

    # flxs
    plot_res(;
        (;
            dat_prefix = :flx,
            xlabel = "model flux",
            ylabel = "exp flux",
            norm_lims = [-1.0, 1.0],
            lims = [D[:mflx] , D[:Mflx]],
            fname = "$(fileid)_tot_corr__flxs.png"
        )...
    )

    plot_res(;
        (;
            dat_prefix = :conc,
            xlabel = "model conc",
            ylabel = "exp conc",
            norm_lims = [0.0, 1.0],
            lims = [D[:mconc], D[:Mconc]],
            fname = "$(fileid)_tot_corr__conc.png"
        )...
    )
end

## -------------------------------------------------------------------
# Error histogram
let
    flx_glc = D[:fba][(:flx, "GLC")]
    # conc_glc = Hd.val("cGLC")

    function _plot(k; d = 0.1)
        p = plot(xlabel = "abs ( exp - model ) / exp ", ylabel = "prob density")
        diff = abs.(D[k][(:flx, :all)] .- D[:Hd][(:flx, :all)]) 
        norm = D[:Hd][(:flx, :all)] 
        ndiff = diff ./ norm
        m, M = extrema(ndiff)
        histogram!(p, diff ./ norm; 
            bins = floor(Int, (M - m)/ d), 
            normalize = :pdf, label = string(k), 
            color = :black, xlim = [0,3]
        )
    end
    p = plot(_plot(:ep), _plot(:fba), layout = 2)
    mysavefig(p, "$(fileid)_model_exp_err_hitograms.png")
end
