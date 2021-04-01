## ----------------------------------------------------------------------------
# beta vs stuff
let
    method = ME_Z_EXPECTED_G_MOVING
    cGLC_plt = plot(;xlabel = "cGLC", ylabel = "beta")
    D_plt = plot(;xlabel = "D", ylabel = "beta")
    for exp in EXPS 
        datfile = ME_INDEX[method, :DFILE, exp]
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

## ----------------------------------------------------------------------------
# flux vs beta
let
    method = ME_Z_EXPECTED_G_MOVING
    p = plot(title = iJR.PROJ_IDER, xlabel = "beta", ylabel = "biom")
    for exp in EXPS 
        datfile = ME_INDEX[method, :DFILE, exp]
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
