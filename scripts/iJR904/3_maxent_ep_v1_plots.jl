using DrWatson 
quickactivate(@__DIR__, "Chemostat_Heerden2013")

using SparseArrays
using Serialization
using StatsBase
using Measures
using MathProgBase
using Statistics

## -------------------------------------------------------------------
# run add "https://github.com/josePereiro/Chemostat" in the 
# julia Pkg REPL for installing the package
import Chemostat
Ch = Chemostat


# Run add https://github.com/josePereiro/Chemostat_Heerden2013.jl in the Julia Pkg REPL to install the
# package, then you must activate the package enviroment (see README)
import Chemostat_Heerden2013: HeerdenData, BegData, iJR904, save_data, load_data
import Chemostat_Heerden2013.iJR904: OBJ_IDER, ATPM_IDER, COST_IDER
const Hd  = HeerdenData;
const iJR = iJR904

## -------------------------------------------------------------------
notebook_name = "2_maxent_ep_v1_plots";

## -------------------------------------------------------------------
# LOADING RESULTS
bundles = load_data(iJR.MAXENT_FBA_EB_BOUNDLES_FILE)
# cache_file = joinpath(iJR.MODEL_PROCESSED_DATA_DIR, "maxent_ep_v1__cache4.jls")
# params, bundles = deserialize(cache_file)
# # bundles = deserialize(cache_file)
# println(realpath(cache_file), " loaded!!!")
exps = bundles |> keys |> collect |> sort
for (exp, bundle) in bundles
    sort!(bundle.βs)
    sort!(bundle.ξs)
end

## -------------------------------------------------------------------
# FIND EXP BETA
# Find the beta that better approximate the objective function with respect 
# to the experimetal data. See Cossio's paper (see README)
closest_βs = Dict()
for exp in eachindex(Hd.val("cGLC"))
    bundle = bundles[exp]
    ξ = Hd.val("xi", exp)
    closest_β = bundle.βs[1]
    for β in bundle.βs
        μ = Ch.Utils.av(bundle, ξ, β, :ep, OBJ_IDER)
        cμ = Ch.Utils.av(bundle, ξ, closest_β, :ep, OBJ_IDER)
        exp_μ = Hd.val("D", exp)
        
        abs(μ - exp_μ) < abs(cμ - exp_μ) && (closest_β = β)
    end
    closest_βs[exp] = closest_β
end

## -------------------------------------------------------------------
# PLOTS
using Plots
# pyplot()

## -------------------------------------------------------------------
color_pool = [:orange, :blue, :red, :black, :violet, 
    :gray, :green, :brown, :magenta];
colors = Dict(exp => rand(color_pool) for exp in exps)
marker_pool = [:circle, :star, :square, :dot]
markers = Dict(exp => rand(marker_pool) for exp in exps)


## -------------------------------------------------------------------
# TODO package this
# CONC
# Transforme a exchange value to a medium concentration (mM), it is by using the 
# steady state assumption in the chemostat, see Cossio's paper (see README)
function conc(Hd_met, exp, ξ::Real, data_idxs...)
    # s = c - u*ξ if u > 0 means intake
    model_met = iJR.Hd_mets_map[Hd_met]
    model_exch = iJR.exch_met_map[model_met]
    bundle = bundles[exp]
    
    u = Ch.Utils.av(bundle, ξ, data_idxs..., model_exch) # exchange val
    uerr = Ch.Utils.va(bundle, ξ, data_idxs..., model_exch) |> sqrt # exchange std
    c = Hd_met == "GLC" ? Hd.val("cGLC", exp) : 0.0
    return (max(c + u * ξ, 0.0), uerr * ξ)
end

function conc(Hd_met, exp, ξs::Vector, data_idxs...)
    model_concs, model_conc_stds = [], []
    foreach(function (ξ) 
                model_conc, model_conc_std = conc(Hd_met, exp, ξ, data_idxs...)
                push!(model_concs, model_conc)
                push!(model_conc_stds, model_conc_std)
            end, ξs)
    return model_concs, model_conc_stds
end

## -------------------------------------------------------------------
# Closest beta
# Just for checking that the experimental objective is inside the beta intervals
# and evaluate the 'experimental' beta approximation
p = Plots.plot(xlabel = "beta", ylabel = "mu", legend = :topleft)
for exp in eachindex(Hd.val("cGLC"))
    
    # model
    exp_ξ = Hd.val("xi", exp)
    bundle = bundles[exp]
    μs = Ch.Utils.av(bundle, exp_ξ, bundle.βs, :ep, OBJ_IDER)
    c = rand(colors)
    m = rand(markers)
    Plots.plot!(p, bundle.βs, μs, label = "", color = colors[exp], lw = 3)

    # exp
    Plots.scatter!(p, [closest_βs[exp]], [Hd.val("D", exp)], label = "", 
        marker = markers[exp], color = colors[exp], ms = 9)
end
p

## -------------------------------------------------------------------
# MEDIUM CONCENTRATION VS XI
# Metabolites
ps = []
for Hd_ider in Hd.msd_mets
    try
        p = Plots.plot(title = Hd_ider, xlabel = "xi", ylabel = "medium conc (mM)")
        for (exp, bundle) in bundles
            
            # FBA
            concs_, conc_errs_ = conc(Hd_ider, exp, bundle.ξs, :fba)
            Plots.plot!(p, bundle.ξs, concs_, color = colors[exp], label = "", ls = :dash)
            
            # EP
            β = closest_βs[exp]
            concs_, conc_errs_ = conc(Hd_ider, exp, bundle.ξs, β, :ep)
            Plots.plot!(p, bundle.ξs, concs_, 
                yerr = conc_errs_,
                color = colors[exp], label = "", 
                m = :star)
            
            Plots.scatter!(p, [Hd.val("xi", exp)], [Hd.val("s$Hd_ider", exp)], 
                color = colors[exp], label = "", ms = 10)
        end
        Plots.plot!(p, [], [], color = :black, ls = :dash, label = "FBA", lw = 3)
        Plots.plot!(p, [], [], color = :black, label = "EP", lw = 3, m = :star)
        Plots.scatter!(p, [], [], color = :black, label = "Exp")
        push!(ps, p)
    catch err
#             println("exp: $exp, met $Hd_met failed!!!") 
    end
end

Plots.plot(ps..., size = [900, 900])

## -------------------------------------------------------------------
# NORMALIZED ERROR
# +
p = Plots.plot(xlabel = "beta", 
    ylabel = "abs_stoi_err/ abs_mean_flx", 
    legend = :topleft)

for (exp, bundle) in bundles
    
    ξ = Hd.val("xi", exp)

    metnet = bundle[ξ, :net]
    
    M,B = length(metnet.mets), length(bundle.βs)
    max_abs_errs_ = []
    mean_abs_errs_ = []
    std_abs_errs_ = []
    
    
    for β in bundle.βs 
        abs_errs_ = abs.(Ch.Utils.stoi_err(bundle, ξ, β, :ep))
        
        max_abs_err_ = maximum(abs_errs_)
        mean_abs_err_ = mean(abs_errs_)
        std_abs_err_ = std(abs_errs_)
        
        
        mean_abs_flxs_ = mean(abs.(Ch.Utils.av(bundle, ξ, β, :ep)))
        
        push!(max_abs_errs_, max_abs_err_/ mean_abs_flxs_)
        push!(mean_abs_errs_, mean_abs_err_/ mean_abs_flxs_)
        push!(std_abs_errs_, std_abs_err_/ mean_abs_flxs_)
    end
    
    Plots.plot!(p, bundle.βs , max_abs_errs_, alpha = 50, color = colors[exp], 
        lw = 1, label = "")
    Plots.plot!(p, bundle.βs , mean_abs_errs_, yerr = std_abs_errs_, alpha = 50, color = colors[exp], 
        lw = 1, label = "", ls = :dash)
    
end
Plots.plot!(p, [], [], lw = 3, label = "max err", color = :black)
Plots.plot!(p, [], [], lw = 3, label = "mean/std err", ls = :dash, color = :black)
p


# ## -------------------------------------------------------------------
# # Stoi error

# function make_stoi_err_gif()
#     p = Plots.plot(title = "iJR904", 
#         xlabel = "beta", ylabel = "abs stoi err", 
#         xaxis = :log10,
#     )

#     for exp in eachindex(Hd.val("cGLC"))

#         bundle = bundles[exp]
#         ξ = Hd.val("xi", exp)
#         ξstr = round(ξ, digits = 2)

#         β = closest_βs[exp]

#         metnet = bundle[ξ, :net)

#         # errs_ = abs.(Ch.Utils.stoi_err(bundle, ξ, β, :ep))
#         # Plots.scatter!(p, fill(β, length(errs_)), [errs_], color = colors[exp], label = "")

#         # flxs_ = [abs.(Ch.Utils.av(bundle, ξ, β, :ep)) for β in bundle.βs]
#         # Plots.plot!(p, bundle.βs, mean.(flxs_), lw = 1, label = "", color = :blue)
#         # Plots.plot!(p, bundle.βs, minimum.(flxs_), lw = 1, label = "", color = :red)
#     end

#     # Legend
#     # Plots.scatter!(p, [], [], color = :black, lw = 3, label = "abs stoi err")
#     # Plots.plot!(p, [], [], ls = :dash, lw = 3, label = "mean abs flx", color = :blue)
#     # Plots.plot!(p, [], [], ls = :dash, lw = 3, label = "min abs flx", color = :red)


#     # gif 
#     # yulims_ = Ch.Utils.logspace(-2.1, 0.3, 40) |> reverse # y axis upper limit
#     # yulims_ = [fill(maximum(yulims_), 25); yulims_]
#     # yulims_ = [reverse(yulims_); yulims_]
#     # yllim_ = -0.0001 # y axis lower bound
#     # gif = Plots.@gif for i in eachindex(yulims_)
#     #     iter_ = yllim_:((yulims_[i] - yllim_)/8):yulims_[i]
#     #     Plots.plot!(p, yaxis = [yllim_, yulims_[i]], yticks = (iter_, 
#     #             string.(round.(collect(iter_), digits = 3))))
#     # end every 1

#     # # saving gif
#     # gif_file = joinpath(iJR.MODEL_FIGURES_DIR, "$(notebook_name)__stoi_err.gif")
#     # cp("tmp.gif", gif_file, force = true)
#     # rm("tmp.gif", force = true)
#     # println(relpath(gif_file), " created!!!")
#     # return gif
# end

# ##
# # gif = make_stoi_err_gif()
# # -

## -------------------------------------------------------------------
# TOTAL CORRELATIONS

ep_p = Plots.plot(xlabel = "exp (mM)", ylabel = "model (mM)", title = "EP")
fba_p = Plots.plot(xlabel = "exp (mM)", ylabel = "model (mM)", title = "FBA")
for exp in eachindex(Hd.val("cGLC"))
    for Hd_met in Hd.msd_mets
        try
            exp_conc = Hd.val("s$Hd_met", exp)
            
            # EP
            model_conc, model_conc_err = conc(Hd_met, exp, Hd.val(:xi, exp), closest_βs[exp], :ep)
            Plots.scatter!(ep_p, [exp_conc], [model_conc], 
                label = "", color = colors[exp], yerr = model_conc_err, ms = 8)
            
            # FBA
            model_conc, model_conc_err = conc(Hd_met, exp, Hd.val(:xi, exp), :fba)
            Plots.scatter!(fba_p, [exp_conc], [model_conc], marker = :star,
                label = "", color = colors[exp], ls = :dash, ms = 10)
        catch err 
#             println("exp: $exp, met $Hd_met failed!!!") 
        end
    end
end
lrange = [-10, 250]
Plots.plot!(ep_p, lrange, lrange; 
    lw = 3, ls = :dash, color = :black, label = "", alpha = 50)
Plots.plot!(fba_p, lrange, lrange; 
    lw = 3, ls = :dash, color = :black, label = "", alpha = 50)
plot(ep_p, fba_p, size = [600, 300])

## -------------------------------------------------------------------
# INDIVIDUAL CORRELATIONS

# +
# Biomass
ps = []
ider = OBJ_IDER
ep_p = Plots.plot(xlabel = "exp (1/ h)", ylabel = "model (1/ h)", title = "EP ($ider)")
fba_p = Plots.plot(xlabel = "exp (1/ h)", ylabel = "model (1/ h)", title = "FBA ($ider)")
for (exp, bundle) in bundles
    β = closest_βs[exp]
    
    exp_av = Hd.val("D", exp)
    
    # EP
    model_av = Ch.Utils.av(bundle, Hd.val(:xi, exp), closest_βs[exp], :ep, ider)
    model_va = Ch.Utils.va(bundle, Hd.val(:xi, exp), closest_βs[exp], :ep, ider) .|> sqrt
    
    Plots.scatter!(ep_p, [exp_av], [model_av], 
                label = "", color = colors[exp], yerr = model_va, ms = 5)
    
    # FBA
    model_av = Ch.Utils.av(bundle,Hd.val(:xi, exp), :fba, ider)
    Plots.scatter!(fba_p, [exp_av], [model_av], marker = :star,
                label = "", color = colors[exp], ls = :dash, ms = 5)
    
end
# Plots.plot!(ep_p, x -> x, lw = 3, ls = :dash, color = :black, label = "", alpha = 50)
# Plots.plot!(fba_p, x -> x, lw = 3, ls = :dash, color = :black, label = "", alpha = 50)
p = plot(ep_p, fba_p, size = [600, 300])
push!(ps, p);

for Hd_met in Hd.msd_mets
    try
        ep_p = Plots.plot(xlabel = "exp (mM)", ylabel = "model (mM)", title = "EP ($Hd_met)")
        fba_p = Plots.plot(xlabel = "exp (mM)", ylabel = "model (mM)", title = "FBA ($Hd_met)")

        for exp in eachindex(Hd.val("cGLC"))
            exp_conc = Hd.val("s$Hd_met", exp)
            
            # EP
            model_conc, model_conc_err = conc(Hd_met, exp, Hd.val(:xi, exp), closest_βs[exp], :ep)
            Plots.scatter!(ep_p, [exp_conc], [model_conc], 
                label = "", color = colors[exp], yerr = model_conc_err, ms = 5)
            
            # FBA
            model_conc, model_conc_err = conc(Hd_met, exp, Hd.val(:xi, exp), :fba)
            Plots.scatter!(fba_p, [exp_conc], [model_conc], marker = :star,
                label = "", color = colors[exp], ls = :dash, ms = 5)
        end
        # Plots.plot!(ep_p, x -> x, lw = 3, ls = :dash, color = :black, label = "", alpha = 50)
        # Plots.plot!(fba_p, x -> x, lw = 3, ls = :dash, color = :black, label = "", alpha = 50)
        p = plot(ep_p, fba_p, size = [600, 300])
        push!(ps, p)

    catch err 
#             println("exp: $exp, met $Hd_met failed!!!") 
    end
end

p = Plots.plot(ps..., size = [800,1600], 
    layout = grid(length(ps), 1), margin = 8mm, titlefont = 13, guidefont = 13)

# Plots.savefig(p, "correlations.pdf")

## -------------------------------------------------------------------
# # ### flx vs xi

# # Biomass
# ider = params["obj_ider"]
# p = Plots.plot(xlabel = "xi", ylabel = "flx", title = ider)
# for (exp, bundle) in bundles
#     β = closest_βs[exp]
    
#     ep_avs = Ch.Utils.av(bundle, bundle.ξs, β, :ep, ider)
#     ep_stds = Ch.Utils.va(bundle, bundle.ξs, β, :ep, ider) .|> sqrt
    
#     fba_avs = Ch.Utils.av(bundle, bundle.ξs, :fba, ider)
    
#     Plots.plot!(p, bundle.ξs, ep_avs, color = colors[exp], label = "")
#     Plots.plot!(p, bundle.ξs, fba_avs, ls = :dash, color = colors[exp], label = "")
#     Plots.scatter!(p, [Hd.val("xi", exp)], [Hd.val("D", exp)], color = colors[exp], label = "")
# end
# p

# # tot_cos
# ider = "tot_cost"
# p = Plots.plot(xlabel = "xi", ylabel = "flx", title = ider)
# for (exp, bundle) in bundles
#     β = closest_βs[exp]
    
#     ep_avs = Ch.Utils.av(bundle, bundle.ξs, β, :ep, ider)
#     ep_stds = Ch.Utils.va(bundle, bundle.ξs, β, :ep, ider) .|> sqrt
    
#     fba_avs = Ch.Utils.av(bundle, bundle.ξs, :fba, ider)
    
#     Plots.plot!(p, bundle.ξs, ep_avs, color = colors[exp], label = "")
#     Plots.plot!(p, bundle.ξs, fba_avs, ls = :dash, color = colors[exp], label = "")
#     Plots.scatter!(p, [Hd.val("xi", exp)], [Hd.val("D", exp)], color = colors[exp], label = "")
# end
# p

# # Measured metabolites
# Hd_met = Hd.msd_mets[3]
# p = Plots.plot(xlabel = "xi", ylabel = "conc", title = Hd_met, xaxis = [0.0, 100])
# for (exp, bundle) in bundles
#     β = closest_βs[exp]
    
#     # ep
#     model_concs, model_conc_stds = conc(Hd_met, exp, bundle.ξs, β, :ep)
#     Plots.plot!(p, bundle.ξs, model_concs, yerr = model_conc_stds, color = colors[exp], label = "")
    
#     # fba
#     model_concs, model_conc_stds = conc(Hd_met, exp, bundle.ξs, :fba)
#     Plots.plot!(p, bundle.ξs, model_concs, ls = :dash, color = colors[exp], label = "")
    
#     Plots.scatter!(p, [Hd.val("xi", exp)], [Hd.val("s$Hd_met", exp)], color = colors[exp], label = "", ms =8)
# end
# p

# # ### Marginals

# ## Biomass
# ider = params["obj_ider"]
# p = Plots.plot(xlabel = "flx", ylabel = "pdf", title = ider, 
#     xaxis = [-0.02, 0.32])
# for (exp, bundle) in bundles |> sort
#     ξ = Hd.val("xi", exp)
#     β = closest_βs[exp]
#     Ch.Plots.plot_marginal!(p, bundle, ξ, β, [:ep], ider, color = colors[exp], label = false)
#     vline!([Hd.val("D", exp)], label = "", ls = :dash, color = colors[exp], lw = 2)
#     Plots.plot!(p, [], [], label = exp, lw = 5, legendtitle = "exps")
# end
# p

# ## rand reaction
# model = deserialize(iJR.BASE_MODEL_FILE)
# ider = rand(model.rxns)
# p = Plots.plot(xlabel = "flx", ylabel = "pdf", title = ider)
# for (exp, bundle) in bundles |> sort
#     ξ = Hd.val("xi", exp)
#     β = closest_βs[exp]
#     Ch.Plots.plot_marginal!(p, bundle, ξ, β, [:ep, :fba], ider,  label = false)
#     Plots.plot!(p, [], [], label = exp, lw = 5, legendtitle = "exps")
# end
# p

# # ### Raw data

# ps = []
# for id in Hd.data_ids()
#     id == "cGLC" && continue
    
#     xs = Hd.val("cGLC")
#     xunit_ = Hd.unit("cGLC")
#     xmin_, xmax_ = minimum(xs), maximum(xs)
#     margin_ = (xmax_ - xmin_)/ 10
#     xlim_ = [xmin_ - margin_, xmax_ + margin_]
    
#     ys = Hd.val(id)
#     yunit_ = Hd.unit(id)
#     ymin_, ymax_ = minimum(ys), maximum(ys)
#     margin_ = (ymax_ - ymin_)/ 10
#     ylim_ = [ymin_ - margin_, ymax_ + margin_]
    
#     p = Plots.scatter(xs, ys, 
#         xlabel = "cGLC ($xunit_)", xaxis = xlim_,
#         ylabel = "$id ($yunit_)", yaxis = ylim_,
#         label = "", title = "")
#     push!(ps, p)
# end

# Plots.plot(ps..., size = [1000, 900], font = 8)


