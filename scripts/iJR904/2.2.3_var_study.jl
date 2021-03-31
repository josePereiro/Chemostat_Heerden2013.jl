## -------------------------------------------------------------------
# Var vs sg
let
    method = ME_MAX_POL

    exp = 5
    datfile = INDEX[method, :DFILE, exp]
    dat = deserialize(datfile)
    
    objider = iJR.BIOMASS_IDER
    model = dat[:model]
    objidx = ChU.rxnindex(model, objider)
    exp_beta = dat[:exp_beta] 
    epout = dat[:epouts][exp_beta]
    exp_xi = Hd.val(:xi, exp)

    vas = log2.(ChU.va(epout))

    p = scatter(;title = "variance histogram", xlabel = "log var")
    histogram!(p, vas)
    mysavefig(p, "ep_variance_hist")
end
