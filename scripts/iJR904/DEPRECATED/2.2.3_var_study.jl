let
    method = ME_MAX_POL
    niders = length(FLX_IDERS)

    stds = map(EXPS) do exp
        DAT[method, :eperr, :flx, FLX_IDERS, exp] ./ abs.(DAT[method, :Hd, :flx, "GLC", exp])
    end
    
    # sGLC
    ids = [:D, :sAC, :sGLC, :sSUCC, :sFORM,  :DCW, :uAC, :uGLC, :uSUCC, :uFORM]
    id_valss = Dict(id => [Hd.val(id, exp) for exp in EXPS] for id in ids)
    push!(ids, :cgD_X)
    id_valss[:cgD_X] = [-Hd.ciD_X(:GLC, exp) for exp in EXPS]

    for id in ids
        p = scatter(;xlabel = string(id), ylabel = "std")
        id_vals = id_valss[id]
        sidx = sortperm(id_vals)
        for idx in sidx
            x = id_vals[idx]
            ys = stds[idx]
            scatter!(p, fill(x, niders), ys; label = "", color = :black, m = 8)
            scatter!(p, [x], [mean(ys)]; label = "", color = :red, m = 12)
        end
        mysavefig(p, "av_var_vs_"; id)
    end
end    