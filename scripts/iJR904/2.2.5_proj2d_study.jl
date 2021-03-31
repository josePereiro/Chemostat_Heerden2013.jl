
## ----------------------------------------------------------------------------
# proj 2D
let
    method = ME_MAX_POL
    biom_ider = iJR.BIOMASS_IDER

    ps_pool = Dict()
    for exp in EXPS

        datfile = INDEX[method, :DFILE, exp]
        dat = deserialize(datfile)
        
        model = dat[:model]
        
        for Hd_ider in FLX_IDERS

            # 2D Projection
            p = plot(;title = string("Heerden2013, exp", exp), 
                xlabel = string(biom_ider), ylabel = string(Hd_ider),
                legend = :right
            )
            proj = DAT[method, :ep, :proj, Hd_ider, exp]
            ChP.plot_projection2D!(p, proj; l = 50)

            # cgD/X
            input = -Hd.cval(Hd_ider, exp, 0.0) * Hd.val(:D, exp) / Hd.val(:DCW, exp)
            hline!(p, [input]; lw = 3, color = :black, ls = :solid, label = "input")

            # EXPERIMENTAL FLXS
            exp_biom = DAT[method, :Hd, :flx, biom_ider, exp]
            exp_exch = DAT[method, :Hd, :flx, Hd_ider, exp]
            scatter!(p, [exp_biom], [exp_exch]; 
                m = 8, color = :red, label = "exp"
            )
            
            # MAXENT FLXS
            ep_biom = DAT[method, :ep, :flx, biom_ider, exp]
            ep_biom_err = DAT[method, :eperr, :flx, biom_ider, exp]
            ep_exch = DAT[method, :ep, :flx, Hd_ider, exp]
            ep_exch_err = DAT[method, :eperr, :flx, Hd_ider, exp]
            scatter!(p, [ep_biom], [ep_exch]; 
                xerr = [ep_biom_err], yerr = [ep_exch_err],
                m = 8, color = :blue, label = "maxent"
            )

            # mysavefig(p, "polytope"; Hd_ider, exp, method)
            ps_pool[(exp, Hd_ider)] = deepcopy(p)
        end
    end

    # collect 
    for exp in EXPS
        ps = Plots.Plot[ps_pool[(exp, Hd_ider)] for Hd_ider in FLX_IDERS]
        mysavefig(ps, "polytope"; exp, method)
    end

    for Hd_ider in FLX_IDERS
        ps = Plots.Plot[ps_pool[(exp, Hd_ider)] for exp in EXPS]
        mysavefig(ps, "polytope"; Hd_ider, method)
    end
end