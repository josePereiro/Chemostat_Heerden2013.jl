
## ----------------------------------------------------------------------------
# marginal distributions
let 
    objider = iJR.BIOMASS_IDER
    size = [300, 250]
    method2 = ME_MAX_POL

    # Iders
    model_iders, Hd_iders = [iJR.BIOMASS_IDER], ["D"]
    for Hd_met in CONC_IDERS
        model_met = Hd_mets_map[Hd_met]
        model_exch = Hd_rxns_map[Hd_met]
        push!(model_iders, model_exch)
        push!(Hd_iders, string("u", Hd_met))
    end
    
    for (model_ider, Hd_ider) in zip(model_iders, Hd_iders)
        ps = Plots.Plot[]
        ps2 = Plots.Plot[]
        for exp in EXPS
            p = plot(title = string(Hd_ider, " exp: ", exp))
            p2 = plot(title = string(Hd_ider, " exp: ", exp))
            margin, m, M = -Inf, Inf, -Inf
            Hd_av = Hd.val(Hd_ider, exp)
            
            # EP
            for method in ALL_METHODS
                color = method_colors[method]    

                datfile = ME_INDEX[method, :DFILE, exp]
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

                if method == method2
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
            xticks =  (EXPS, string.(EXPS))
            p = bar!(p, EXPS, Hd.val(k); title = k, label = "", xticks)
            push!(ps, p)
            push!(ps2, p)
        end

        pname = string(Hd_ider, "_marginals")
        mysavefig(ps, pname)
        pname = string(Hd_ider, "_marginals_vs_beta")
        mysavefig(ps2, pname)
    end

end 

## ----------------------------------------------------------------------------
# marginals v2
let 
    objider = iJR.BIOMASS_IDER
    size = [300, 250]

    # Iders
    model_iders, Hd_iders = [iJR.BIOMASS_IDER], ["D"]
    for Hd_met in CONC_IDERS
        model_met = Hd_mets_map[Hd_met]
        model_exch = Hd_rxns_map[Hd_met]
        push!(model_iders, model_exch)
        push!(Hd_iders, string("u", Hd_met))
    end
    
    for (model_ider, Hd_ider) in zip(model_iders, Hd_iders)
        marg_params = (;xlabel = string(Hd_ider), yaxis = nothing, ylabel = "prob", titlefont = 8)

        epps = Plots.Plot[]
        exps = Plots.Plot[]
        for method in ALL_METHODS
            expp = plot(;title = string("Experimental"), marg_params...)
            epp = plot(;title = string(" MaxEnt: \n", method), marg_params...)
            margin, m, M = -Inf, Inf, -Inf
            
            # EP
            for exp in EXPS
                Hd_av = Hd.val(Hd_ider, exp)
                color = exp_colors[exp]    

                datfile = ME_INDEX[method, :DFILE, exp]
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
                xlabel = "rep", ylabel = string(k), titlefont = 8
            )
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
