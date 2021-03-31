## ----------------------------------------------------------------------------
# EP biomass corr
let
    ps = Plots.Plot[]
    for method in ALL_METHODS
        p = plot(title = string(iJR.PROJ_IDER, " method: ", method), 
            xlabel = "exp biom", ylabel = "model biom")
        ep_vals = DAT[method, :ep, :flx, iJR.BIOMASS_IDER, EXPS]
        eperr_vals = DAT[method, :eperr, :flx, iJR.BIOMASS_IDER, EXPS]
        Hd_vals = DAT[method, :Hd, :flx, iJR.BIOMASS_IDER, EXPS]
        m, M = myminmax([Hd_vals; ep_vals])
        margin = abs(M - m) * 0.1
        scatter!(p, Hd_vals, ep_vals; 
            yerr = eperr_vals,
            label = "", color = :black,
            alpha = 0.7, ms = 7,
            xlim = [m - margin, M + margin],
            ylim = [m - margin, M + margin],
        )
        push!(ps, p)
    end
    layout = (1, length(ps))
    mysavefig(ps, "obj_val_ep_corr"; layout)
end

## ----------------------------------------------------------------------------
# total correlations
let
    for (dat_prefix, iders) in [(:flx, FLX_IDERS), (:conc, CONC_IDERS)]

        ps = Plots.Plot[]
        for method in ALL_METHODS                                            
            ep_vals = DAT[method, :ep, dat_prefix, iders, EXPS]
            ep_errs = DAT[method, :eperr, dat_prefix, iders, EXPS]
            Hd_vals = DAT[method, :Hd, dat_prefix, iders, EXPS]
            color = [ider_colors[ider] for ider in iders, exp in EXPS]
            
            diffsign = sign.(Hd_vals) .* sign.(ep_vals)
            Hd_vals = abs.(Hd_vals) .* diffsign
            ep_vals = abs.(ep_vals) .* diffsign
            
            m, M = myminmax([ep_vals; Hd_vals])
            scatter_params = (;label = "", color, ms = 7, alpha = 0.7)
            # ep corr
            p1 = plot(title = "$(iJR.PROJ_IDER) (EP) $method", 
                ylabel = "model signdiff $(dat_prefix)",
                xlabel = "exp signdiff $(dat_prefix)",
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
