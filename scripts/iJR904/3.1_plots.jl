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

using Plots

## -------------------------------------------------------------------
# save results
DATA = ChU.load_data(iJR.MAXENT_FBA_EB_BOUNDLES_FILE)
SDATA = sort(collect(DATA); by = first);

## -------------------------------------------------------------------
# flux vs beta
let
    p = plot(title = iJR.PROJ_IDER, xlabel = "beta", ylabel = "biom")
    for (exp, D) in SDATA
        model = D["model"]
        objidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
        epouts = D["epouts"]
        exp_beta = D["exp_beta"]
        exp_xi = D["exp_xi"]
        scatter!(p, [exp_beta], [Hd.val("D", exp)], ms = 12, color = :white, label = "")

        betas = collect(keys(epouts)) |> sort
        bioms = [ChU.av(model, epouts[beta], objidx) for beta in betas]
        scatter!(p, betas, bioms, label = "", color = :black, alpha = 0.2)

    end
    p
end

## -------------------------------------------------------------------
# beta_exp corr
let
    p = plot(title = iJR.PROJ_IDER, xlabel = "exp biom", ylabel = "model biom")
    for (exp, D) in SDATA

        model = D["model"]
        objidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
        exp_beta = D["exp_beta"]
        epout = D["epouts"][exp_beta]

        scatter!(p, [ChU.av(model, epout, objidx)], [Hd.val("D", exp)], 
            color = :white, label = "")

    end
    p
end

## -------------------------------------------------------------------
# total correlations
let
    Hd_vals = Dict()
    fba_vals = Dict()
    ep_vals = Dict()
    ep_stds = Dict()
    objider = iJR.BIOMASS_IDER
    mflx, Mflx = Inf, -Inf
    mconc, Mconc = Inf, -Inf

    getpush!(d, kvs...) = foreach((kv) -> push!(get!(d, kv[1], []), kv[2]), kvs)

    for (exp, D) in SDATA

        model = D["model"]
        objidx = ChU.rxnindex(model, objider)
        exp_beta = D["exp_beta"]
        epout = D["epouts"][exp_beta]
        exp_xi = Hd.val("xi", exp)
        fbaout = D["fbaout"]

        # Biomass
        fba_biom = ChU.av(model, fbaout, objidx)
        ep_biom = ChU.av(model, epout, objidx)
        ep_std = sqrt(ChU.va(model, epout, objidx))
        Hd_biom = Hd.val("D", exp)
        
        getpush!(fba_vals, ((:flx, objider), fba_biom))
        getpush!(ep_vals, ((:flx, objider), ep_biom))
        getpush!(ep_stds, ((:flx, objider), ep_std))
        getpush!(Hd_vals, ((:flx, objider), Hd_biom))
        
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
                flxkey = (:flx, Hd_met)
                conckey = (:conc, Hd_met)
                getpush!(fba_vals, (flxkey, fba_av), (conckey, fba_conc))
                getpush!(ep_vals, (flxkey, ep_av), (conckey, ep_conc))
                getpush!(ep_stds, (flxkey, ep_std), (conckey, ep_std * exp_xi))
                getpush!(Hd_vals, (flxkey, Hd_flx), (conckey, Hd_conc))

                mflx = minimum([mflx, fba_av, ep_av, fba_av])
                Mflx = maximum([Mflx, fba_av, ep_av, fba_av])
                mconc = minimum([mconc, fba_conc, ep_conc, Hd_conc])
                Mconc = maximum([Mconc, fba_conc, ep_conc, Hd_conc])
                
            catch 
                @warn string(Hd_met, " fails")
            end
        end

    end

    # plots
    function plot_res(;xlabel, ylabel, dat_prefix, norm_lims, lims, fname)

        ep_p = plot(;title = "$(iJR.PROJ_IDER) (EP)", xlabel, ylabel)
        ep_p_norm = plot(;title = "$(iJR.PROJ_IDER) (EP normalized)", xlabel, ylabel)
        fba_p = plot(;title = "$(iJR.PROJ_IDER) (FBA)", xlabel, ylabel)
        fba_p_norm = plot(;title = "$(iJR.PROJ_IDER) (FBA normalized)", xlabel, ylabel)

        for ider in [Hd.msd_mets; objider]

            # flx
            flxkey = (dat_prefix, ider)
            if haskey(Hd_vals, flxkey)
            
                scatter!(fba_p, Hd_vals[flxkey], fba_vals[flxkey]; color = :black, label = "")
                norm = maximum([abs.(fba_vals[flxkey]); abs.(Hd_vals[flxkey])])
                scatter!(fba_p_norm, Hd_vals[flxkey] ./ norm, fba_vals[flxkey] ./ norm; color = :black, label = "")
                
                scatter!(ep_p, Hd_vals[flxkey], ep_vals[flxkey];
                    yerr = ep_stds[flxkey], color = :black, label = "")
                norm = maximum([abs.(ep_vals[flxkey]); abs.(Hd_vals[flxkey])])
                scatter!(ep_p_norm, Hd_vals[flxkey] ./ norm, ep_vals[flxkey] ./ norm; 
                    yerr = ep_stds[flxkey] ./ norm, color = :black, label = "")
            end
        end
        plot!(fba_p_norm, norm_lims, norm_lims; ls = :dash, color = :black, label = "")
        plot!(ep_p_norm, norm_lims, norm_lims; ls = :dash, color = :black, label = "")
        
        plot!(fba_p, lims, lims; ls = :dash, color = :black, label = "")
        plot!(ep_p, lims, lims; ls = :dash, color = :black, label = "")
        

        p = plot([ep_p_norm, fba_p_norm, ep_p, fba_p]...; layout = 4, size = [800, 800])
        fig_path = joinpath(iJR.MODEL_FIGURES_DIR, fname)
        savefig(p, fig_path)
        
        @info "Saved at" fig_path
    end

    # flxs
    plot_res(;
        dat_prefix = :flx,
        xlabel = "model flux",
        ylabel = "exp flux",
        norm_lims = [-1.0, 1.0],
        lims = [mflx , Mflx],
        fname = "tot_corr__flxs.png"
    )

    plot_res(;
        dat_prefix = :conc,
        xlabel = "model conc",
        ylabel = "exp conc",
        norm_lims = [0.0, 1.0],
        lims = [mconc, Mconc],
        fname = "tot_corr__conc.png"
    )
end
