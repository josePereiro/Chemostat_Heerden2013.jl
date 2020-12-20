import DrWatson: quickactivate, savename
quickactivate(@__DIR__, "Chemostat_Heerden2013")

import SparseArrays

# -------------------------------------------------------------------
# Run add https://github.com/josePereiro/Chemostat_Heerden2013.jl in the Julia Pkg REPL to install the
# package, then you must activate the package enviroment (see README)
import Chemostat_Heerden2013
const Hd  = Chemostat_Heerden2013.HeerdenData;
const BD  = Chemostat_Heerden2013.BegData;
const iJR = Chemostat_Heerden2013.iJR904

# -------------------------------------------------------------------
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
using Plots.Measures
using Serialization
using Statistics

## -------------------------------------------------------------------
var_order = 49
fname = savename("exploration_", (;var_order), "jld")
fpath = joinpath(iJR.MODEL_PROCESSED_DATA_DIR, fname)
D = deserialize(fpath)
@info "Data info" var_order length(D["epouts"]) D["error"]

model = D["model"]
M, N = size(model)
obj_ider = iJR.BIOMASS_IDER
obj_idx = ChU.rxnindex(model, obj_ider)
glc_ider = "EX_glc_LPAREN_e_RPAREN_"
glc_idx = ChU.rxnindex(model, glc_ider)
epouts = D["epouts"]
sepouts = sort(collect(epouts), by = first)
betas = first.(sepouts)
fbaout = ChLP.fba(model, obj_idx)
ex_glc = fbaout.v[glc_idx]
@assert ex_glc != 0.0
EPS = 1e-8

## -------------------------------------------------------------------
# abs(fba - ep) / ex_glc
let
    p = plot(title = iJR.PROJ_IDER, 
        xlabel = "log10( beta )", 
        ylabel = "log10( abs( fba - ep ) / ex_glc )"
    )
    f(x) = log10(x)
    mat = zeros(length(sepouts), N)
    for (i, (beta, epout)) in enumerate(sepouts)
        mat[i, :] = (abs.(fbaout.v .- epout.av) ./ abs(ex_glc)) .+ EPS
    end
    plot!(p, first.(sepouts), f.(mat); 
        label = "", color = :black, alpha = 0.1
    )
    plot!(p, first.(sepouts), f.(mat[:, obj_idx] .+ EPS);
        label = "", color = :blue, ls = :dash, alpha = 0.8, lw = 3
    )
    p
end

## -------------------------------------------------------------------
# abs(fba - ep) / max_diff
let
    p1 = plot(title = iJR.PROJ_IDER, 
        xlabel = "beta", 
        ylabel = "abs( fba - ep ) / max_diff"
    )
    f(x) = x
    mat = zeros(length(sepouts), N)
    for (i, (beta, epout)) in enumerate(sepouts)
        mat[i, :] = (abs.(fbaout.v .- epout.av) ./ abs(ex_glc)) .+ EPS
    end

    norm = ones(length(sepouts)) * maximum.(eachcol(mat))'
    mat ./= norm

    plot!(p1, first.(sepouts), f.(mat); 
        label = "", color = :black, alpha = 0.1
    )
    plot!(p1, first.(sepouts), f.(mat[:, obj_idx] .+ EPS);
         label = "", color = :blue, ls = :dash, alpha = 0.8, lw = 3
    )

    # histograms
    kwargs = (;
        normalize = :pdf,
        bins = 200,
        xlabel = "prob. density",
        ylabel = "abs( fba - ep ) / max_diff",
        legend = :best, 
        orientation = :h, 
        color = :black
    )
    p2 = histogram(mat[begin, :];
        label = string("beta: ", round(Int, betas[begin])), 
        kwargs...
    )
    p3 = histogram(mat[end, :];  
        label = string("beta: ", round(Int, betas[end])), 
        kwargs...
    )
    plot([p1, p2, p3]..., layout = grid(1, 3), size = [1200, 400], 
        margin = 5mm
    )
end

## -------------------------------------------------------------------
# stoi err
let
    p = plot(title = iJR.PROJ_IDER, 
        xlabel = "log10( beta )", ylabel = "log10( stoi_err / ex_glc )")
    f(x) = log10(x)
    mat = zeros(length(sepouts), M)
    for (i, (beta, epout)) in enumerate(sepouts)
        mat[i, :] = (abs.(ChU.stoi_err(model, epout)) ./ abs(ex_glc)) .+ EPS
    end
    plot!(p, first.(sepouts), f.(mat); 
        label = "", color = :black, alpha = 0.1
    )
    plot!(p, first.(sepouts), f.(mean.(eachrow(mat))); 
        label = "mean", color = :red, alpha = 0.8, 
        lw = 5, ls = :dash
    )
    p
end