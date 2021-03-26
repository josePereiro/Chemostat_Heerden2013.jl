import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

@time begin
    
    import Chemostat_Heerden2013
    const ChH = Chemostat_Heerden2013
    const Hd = ChH.HeerdenData    

    using UtilsJL
    const UJL = UtilsJL

    using Plots
end

## -------------------------------------------------------------------
fileid = "1.1"
function mysavefig(p, pname; params...) 
    fname = UJL.mysavefig(p, string(fileid, "_", pname), Hd.HEERDEN_FIGURES_DIR; params...)
    @info "Plotting" fname
end

## -------------------------------------------------------------------
# D vs X
let
    p = plot(; title = "Heerden", 
        xlabel = string("Xv (", Hd.unit(:DCW), ")"), 
        ylabel = string("D (", Hd.unit(:D), ")")
    )
    scatter!(p, Hd.val(:D), Hd.val(:DCW); 
        label = "", m = 8, color = :black
    )
    mysavefig(p, "X_vs_D") 
end

## -------------------------------------------------------------------
# BALANCE
let
    ps = Plots.Plot[]
    for met in [:GLC]
        p = plot(; title = string("Balance: ", met), 
            xlabel = "feed", 
            ylabel = "exch + drain" 
        )

        exps = 1:13 

        feed = Hd.cval.(met, exps) .* Hd.val.(:D, exps)
        exch = Hd.uval.(met, exps) .* Hd.val.(:DCW, exps) .|> abs
        drain = Hd.sval.(met, exps) .* Hd.val.(:D, exps) .|> abs

        scatter!(p, feed, exch .+ drain; 
            label = "", m = 8, color = :black
        )
        vals = [feed; exch .+ drain] |> sort
        plot!(p, vals, vals;
            label = "", ls = :dash, lw = 3, alpha = 0.6, color = :black
        )
        push!(ps, p)
    end
    mysavefig(ps, "Balances") 
end