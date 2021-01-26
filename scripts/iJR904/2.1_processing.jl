import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid

    # -------------------------------------------------------------------
    # Run add https://github.com/josePereiro/Chemostat_Heerden2013.jl in 
    # the Julia Pkg REPL to install the package, then you must activate 
    # the package enviroment (see README)
    import Chemostat_Heerden2013
    const ChHd = Chemostat_Heerden2013
    const Hd  = ChHd.HeerdenData;
    const BD  = ChHd.BegData;
    const iJR = ChHd.iJR904

    # -------------------------------------------------------------------
    # run add "https://github.com/josePereiro/Chemostat" in the 
    # julia Pkg REPL for installing the package
    import Chemostat
    import Chemostat.LP.MathProgBase

    const Ch = Chemostat
    const ChP = Ch.Plots
    const ChU = Ch.Utils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP
    const ChEP = Ch.MaxEntEP
    const ChSU = Ch.SimulationUtils

    import UtilsJL
    const UJL = UtilsJL

    using Serialization

    # -------------------------------------------------------------------
    using Plots, FileIO
    import GR
    GR.inline("png")

end

## -------------------------------------------------------------------
INDEX = ChU.load_data(iJR.MAXENT_VARIANTS_INDEX_FILE; verbose = false);

## -------------------------------------------------------------------
const HOMO = :HOMO
const BOUNDED = :BOUNDED
const EXPECTED = :EXPECTED

## -------------------------------------------------------------------
fileid = "2.1"
function mysavefig(p, pname; params...) 
    fname = UJL.mysavefig(p, string(fileid, "_", pname), iJR.MODEL_FIGURES_DIR; params...)
    @info "Plotting" fname
end
myminmax(a::Vector) = (minimum(a), maximum(a))
CONC_IDERS = ["GLC", "SA", "AcA", "FA"]
FLX_IDERS = ["GLC", "SA", "AcA", "FA"]

EXPS = Hd.EXPS # experiments that have both concentration and flx data

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
    HOMO => :red,
    BOUNDED => :orange,
    EXPECTED => :blue,
)

## -------------------------------------------------------------------
DAT = ChU.DictTree()
let 
    objider = iJR.BIOMASS_IDER

    for method in [HOMO, EXPECTED, BOUNDED]
        for exp in Hd.EXPS
            
            datfile = INDEX[method, :DFILE, exp]
            dat = deserialize(datfile)
            
            model = dat[:model]
            objidx = ChU.rxnindex(model, objider)
            exp_beta = dat[:exp_beta]
            epout = dat[:epouts][exp_beta]
            exp_xi = Hd.val(:xi, exp)

            println()
            @info("Doing", exp, method, length(dat[:epouts]), epout.iter); 
            # fbaout = dat[:fbaout]

            # Biomass
            ep_biom = ChU.av(model, epout, objidx)
            ep_std = sqrt(ChU.va(model, epout, objidx))
            Hd_biom = Hd.val("D", exp)
            
            # store
            DAT[method, :ep   , :flx, objider, exp] = ep_biom
            DAT[method, :eperr, :flx, objider, exp] = ep_std
            DAT[method, :Hd   , :flx, objider, exp] = Hd_biom
            DAT[:Hd   , :flx, objider, exp] = Hd_biom
            DAT[method, :fva  , :flx, objider, exp] = ChU.bounds(model, objider)
            
            # mets
            for Hd_met in FLX_IDERS

                    model_met = iJR.Hd_mets_map[Hd_met]
                    model_exch = iJR.exch_met_map[model_met]
                    model_exchi = ChU.rxnindex(model, model_exch)

                    # fuxes
                    ep_av = ChU.av(model, epout, model_exchi)
                    ep_std = sqrt(ChU.va(model, epout, model_exchi))
                    Hd_flx = Hd.val("u$Hd_met", exp)
                    
                    DAT[method, :Hd, :flx, Hd_met, exp] = Hd_flx
                    DAT[:Hd, :flx, Hd_met, exp] = Hd_flx
                    DAT[method, :ep, :flx, Hd_met, exp] = ep_av
                    DAT[method, :eperr, :flx, Hd_met, exp] = ep_std
                    
                    DAT[method, :fva , :flx, Hd_met, exp] = ChU.bounds(model, model_exch)

            end

            for Hd_met in CONC_IDERS

                ep_av = DAT[method, :ep, :flx, Hd_met, exp]
                ep_std = DAT[method, :eperr, :flx, Hd_met, exp] 

                # conc (s = c + u*xi)
                c = Hd.val("c$Hd_met", exp, 0.0)
                ep_conc = max(c + ep_av * exp_xi, 0.0)
                Hd_conc = Hd.val("s$Hd_met", exp)
                
                DAT[method, :Hd, :conc, Hd_met, exp] = Hd_conc
                DAT[:Hd, :conc, Hd_met, exp] = Hd_conc
                DAT[method, :ep, :conc, Hd_met, exp] = ep_conc
                DAT[method, :eperr, :conc, Hd_met, exp] = ep_std * exp_xi
            end

        end # for exp in Hd.EXPS
    end

end

## -------------------------------------------------------------------
# beta vs stuff
let
    method = EXPECTED
    cGLC_plt = plot(;xlabel = "cGLC", ylabel = "beta")
    D_plt = plot(;xlabel = "D", ylabel = "beta")
    for exp in Hd.EXPS 
        datfile = INDEX[method, :DFILE, exp]
        dat = deserialize(datfile)
        beta = dat[:exp_beta]

        params = (;label = "", color = exp_colors[exp], 
            alpha = 0.7, ms = 7
        )
        cGLC = Hd.val("cGLC", exp)
        D = Hd.val("D", exp)
        scatter!(cGLC_plt, [cGLC], [beta]; params...)
        scatter!(D_plt, [D], [beta]; params...)
    end
    mysavefig([cGLC_plt, D_plt], "beta_vs_stuff"; method)
end

## -------------------------------------------------------------------
# EP biomass corr
let
    ps = Plots.Plot[]
    for method in [HOMO, EXPECTED, BOUNDED]
        p = plot(title = string(iJR.PROJ_IDER, " method: ", method), 
            xlabel = "exp biom", ylabel = "model biom")
        ep_vals = DAT[method, :ep, :flx, iJR.BIOMASS_IDER, Hd.EXPS]
        eperr_vals = DAT[method, :eperr, :flx, iJR.BIOMASS_IDER, Hd.EXPS]
        Hd_vals = DAT[method, :Hd, :flx, iJR.BIOMASS_IDER, Hd.EXPS]
        color = [exp_colors[exp] for exp in Hd.EXPS]
        m, M = myminmax([Hd_vals; ep_vals])
        margin = abs(M - m) * 0.1
        scatter!(p, Hd_vals, ep_vals; 
            yerr = eperr_vals,
            label = "", color,
            alpha = 0.7, ms = 7,
            xlim = [m - margin, M + margin],
            ylim = [m - margin, M + margin],
        )
        push!(ps, p)
    end
    layout = (1, length(ps))
    mysavefig(ps, "obj_val_ep_corr"; layout)
end

## -------------------------------------------------------------------
# EXPECTED flux vs beta
let
    method = EXPECTED
    p = plot(title = iJR.PROJ_IDER, xlabel = "beta", ylabel = "biom")
    for exp in Hd.EXPS 
        datfile = INDEX[method, :DFILE, exp]
        dat = deserialize(datfile)
        model = dat[:model]
        objidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
        epouts = dat[:epouts]
        exp_beta = dat[:exp_beta]
        exp_xi = Hd.val("xi", exp)
        scatter!(p, [exp_beta], [Hd.val("D", exp)], ms = 12, color = :white, label = "")

        betas = collect(keys(epouts)) |> sort
        bioms = [ChU.av(model, epouts[beta], objidx) for beta in betas]
        scatter!(p, betas, bioms, label = "", color = :black, alpha = 0.2)

    end
    mysavefig(p, "obj_val_vs_beta"; method)
end

## -------------------------------------------------------------------
# total correlations
let
    for (dat_prefix, iders) in [(:flx, FLX_IDERS), (:conc, CONC_IDERS)]

        ps = Plots.Plot[]
        for method in [HOMO, EXPECTED, BOUNDED]                                            
            ep_vals = DAT[method, :ep, dat_prefix, iders, Hd.EXPS]
            ep_errs = DAT[method, :eperr, dat_prefix, iders, Hd.EXPS]
            Hd_vals = DAT[method, :Hd, dat_prefix, iders, Hd.EXPS]
            color = [ider_colors[ider] for ider in iders, exp in Hd.EXPS]
            
            diffsign = sign.(Hd_vals) .* sign.(ep_vals)
            Hd_vals = abs.(Hd_vals) .* diffsign
            ep_vals = abs.(ep_vals) .* diffsign
            
            m, M = myminmax([ep_vals; Hd_vals])
            scatter_params = (;label = "", color, ms = 7, alpha = 0.7)
            # ep corr
            p1 = plot(title = "$(iJR.PROJ_IDER) (EP) $method", 
                ylabel = "model $(dat_prefix)", 
                xlabel = "exp $(dat_prefix)",
            )
            scatter!(p1, Hd_vals, ep_vals; yerr = ep_errs, scatter_params...)
            plot!(p1, [m,M], [m,M]; ls = :dash, color = :black, label = "")
            push!(ps, deepcopy(p1))

        end

        layout = (1, length(ps))
        pname = string(dat_prefix, "_tot_corr")
        mysavefig(ps, pname; layout)
    end

end

## -------------------------------------------------------------------
# fva bounds
let
   
    ps = Plots.Plot[]
    for ider = FLX_IDERS
        p = plot(title = ider, xlabel = "replica", ylabel = "flx")
        xticks =  (Hd.EXPS, string.(Hd.EXPS))
        
        Hd_vals = DAT[:Hd, :flx, ider, Hd.EXPS]
        plot!(p, Hd.EXPS, Hd_vals; 
            label = "exp", color = :black, alpha = 0.8, lw = 3, xticks)

        for method in [HOMO, EXPECTED, BOUNDED]             
            color = method_colors[method]    
            
            ep_vals = DAT[method, :ep, :flx, ider, Hd.EXPS]
            plot!(p, Hd.EXPS, ep_vals; 
                label = string(method), color, alpha = 0.5, lw = 5, ls = :dash, xticks)
            
            # fba_vals = DAT[method, :fba, :flx, ider, Hd.EXPS]
            # plot!(p, Hd.EXPS, fba_vals; 
            #     label = "fba", color = :red, alpha = 0.5, lw = 5, ls = :dash, xticks)
            
            fva_ranges = DAT[method, :fva, :flx, ider, Hd.EXPS]
            plot!(p, Hd.EXPS, last.(fva_ranges);  
                label = "", color, alpha = 0.8, ls = :dot, lw = 3, xticks)
            plot!(p, Hd.EXPS, first.(fva_ranges); 
                label = "", color, alpha = 0.8, ls = :dot, lw = 3, xticks)
        end
        push!(ps, p)
    end
    pname = string("bound_study")
    mysavefig(ps, pname)
    
end

## -------------------------------------------------------------------
# marginal distributions
let 
    objider = iJR.BIOMASS_IDER
    size = [300, 250]

    # Iders
    model_iders, Hd_iders = [iJR.BIOMASS_IDER], ["D"]
    for Hd_met in CONC_IDERS
        model_met = iJR.Hd_mets_map[Hd_met]
        model_exch = iJR.exch_met_map[model_met]
        push!(model_iders, model_exch)
        push!(Hd_iders, string("u", Hd_met))
    end
    
    for (model_ider, Hd_ider) in zip(model_iders, Hd_iders)
        ps = Plots.Plot[]
        ps2 = Plots.Plot[]
        for exp in Hd.EXPS
            p = plot(title = string(Hd_ider, " exp: ", exp))
            p2 = plot(title = string(Hd_ider, " exp: ", exp))
            margin, m, M = -Inf, Inf, -Inf
            Hd_av = Hd.val(Hd_ider, exp)
            
            # EP
            for method in [BOUNDED, EXPECTED, HOMO]
                color = method_colors[method]    

                datfile = INDEX[method, :DFILE, exp]
                dat = deserialize(datfile)
                model = dat[:model]
                objidx = ChU.rxnindex(model, objider)
                exp_beta = dat[:exp_beta]
                epouts = dat[:epouts]
                epout = epouts[exp_beta]
                ep_av = ChU.av(model, epout, model_ider)
                ep_va = sqrt(ChU.va(model, epout, model_ider))
                # fbaout = dat[:fbaout]
                        
                # ChP.plot_marginal!(p, model, [epout, fbaout], model_exch; legend = false)
                ChP.plot_marginal!(p, model, [epout], model_ider; 
                    legend = false, color, alpha = 0.6, lw = 5)
                
                m = minimum([m, ep_av, Hd_av])
                M = maximum([M, ep_av, Hd_av])
                margin = maximum([margin, 3 * ep_va])

                if method == EXPECTED
                    for (beta, epout) in sort(epouts; by = first)
                        ep_av = ChU.av(model, epout, model_ider)
                        ep_va = sqrt(ChU.va(model, epout, model_ider))
                        
                        alpha = 0.15
                        color = method_colors[method]
                        ChP.plot_marginal!(p2, model, epout, model_ider; 
                            legend = false, color, alpha, lw = 1)

                        if beta == exp_beta
                            # @info "At" Hd_ider exp
                            ChP.plot_marginal!(p2, model, epout, model_ider; 
                                legend = false, color, 
                                alpha = 1.0, lw = 3
                            )
                            break
                        end
                    end
                    push!(ps2, p2)
                end

            end
            # Experimental
            vline!(p, [Hd_av]; label = "", lw = 6, color = :black, alpha = 0.3)
            vline!(p2, [Hd_av]; label = "", lw = 6, color = :black, alpha = 0.3)
            
            plot!(p; xlim = [m - margin, M + margin], size)
            plot!(p2; xlim = [m - margin, M + margin], size)
            push!(ps, p)
        end

        for k in [:xi, :D, :sGLC]
            p = plot(;title = Hd_ider, size)
            xticks =  (Hd.EXPS, string.(Hd.EXPS))
            p = bar!(p, Hd.EXPS, Hd.val(k); title = k, label = "", xticks)
            push!(ps, p)
            push!(ps2, p)
        end

        pname = string(Hd_ider, "_marginals")
        mysavefig(ps, pname)
        pname = string(Hd_ider, "_marginals_vs_beta")
        mysavefig(ps2, pname)
    end

end 

## -------------------------------------------------------------------
# marginals v2
let 
    objider = iJR.BIOMASS_IDER
    size = [300, 250]

    # Iders
    model_iders, Hd_iders = [iJR.BIOMASS_IDER], ["D"]
    for Hd_met in CONC_IDERS
        model_met = iJR.Hd_mets_map[Hd_met]
        model_exch = iJR.exch_met_map[model_met]
        push!(model_iders, model_exch)
        push!(Hd_iders, string("u", Hd_met))
    end
    
    for (model_ider, Hd_ider) in zip(model_iders, Hd_iders)
        marg_params = (;xlabel = string(Hd_ider), yaxis = nothing, ylabel = "prob")

        epps = Plots.Plot[]
        exps = Plots.Plot[]
        for method in [BOUNDED, EXPECTED, HOMO]
            expp = plot(;title = string("Experimental"), marg_params...)
            epp = plot(;title = string(" MaxEnt: ", method), marg_params...)
            margin, m, M = -Inf, Inf, -Inf
            
            # EP
            for exp in EXPS
                Hd_av = Hd.val(Hd_ider, exp)
                color = exp_colors[exp]    

                datfile = INDEX[method, :DFILE, exp]
                dat = deserialize(datfile)
                model = dat[:model]
                objidx = ChU.rxnindex(model, objider)
                epouts = dat[:epouts]
                exp_beta = maximum(keys(epouts))
                epout = epouts[exp_beta]
                ep_av = ChU.av(model, epout, model_ider)
                ep_va = sqrt(ChU.va(model, epout, model_ider))
                
                # bounds = ChU.bounds(model, model_ider)
                ChP.plot_marginal!(epp, model, [epout], model_ider; 
                    legend = false, color, alpha = 0.8, lw = 3
                )
                
                m = minimum([m, ep_av, Hd_av])
                M = maximum([M, ep_av, Hd_av])
                margin = maximum([margin, 3 * ep_va])

                # Experimental
                vline!(expp, [Hd_av]; label = "", lw = 3, color, alpha = 0.8)
                
            end
            
            map([expp, epp]) do p
                plot!(p; xlim = [m - margin, M + margin], size)
            end

            push!(epps, epp)
            push!(exps, expp)
        end

        extras = Plots.Plot[]
        for k in [:xi, :D, :sGLC]
            p = plot(;title = "Experimental", size, 
                xlabel = "rep", ylabel = string(k))
            xticks =  (EXPS, string.(EXPS))
            vals = [Hd.val(k, exp) for exp in EXPS]
            color = [exp_colors[exp] for exp in EXPS]
            p = bar!(p, EXPS, vals; label = "", xticks, color)
            push!(extras, p)
        end

        ps = Plots.Plot[exps; epps; extras]
        layout = (3, 3)
        pname = string(Hd_ider, "_marginals_v2")
        mysavefig(ps, pname; layout)

    end # for (model_ider, Hd_ider)

end 

## -------------------------------------------------------------------
# leyends
# TODO fix this...
let
    for (title, colors) in [
            ("exp", exp_colors), 
            ("iders", ider_colors),
            ("method", method_colors)
        ]
    p = plot(; framestyle = :none)
        scolors = sort(collect(colors); by = (p) -> string(first(p)))
        for (id, color) in scolors
            scatter!(p, [0], [0];
                thickness_scaling = 1,
                color, ms = 8, label = string(id),
                legendfontsize=10, 
                # size = [300, 900],
                # legend = :left
            )
        end
        mysavefig(p, "$(title)_color_legend")
    end
end