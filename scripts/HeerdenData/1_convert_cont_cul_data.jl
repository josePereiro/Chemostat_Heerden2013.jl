import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

## ------------------------------------------------------------------
using CSV
using DataFrames

# Run add https://github.com/josePereiro/Chemostat_Heerden2013.jl in the Julia Pkg REPL to install the
# package, then you must activate the package enviroment (see README)
import Chemostat_Heerden2013: HeerdenData
const Hd = HeerdenData

## ------------------------------------------------------------------
# Load data
# Data taken from Heerden2013 https://doi.org/10.1186/1475-2859-12-80.

# Table 2 Steady state results from continuous fermentations
orig_data = CSV.read(Hd.HEERDEN_CONT_CUL_DATA_ORIG_FILE,  DataFrame;
    delim = "\t", comment = "#");
for n in names(orig_data)
    orig_data[!, n] = Vector{Union{String, Float64}}(orig_data[!, n])
    orig_data[3:end, n] .= parse.(Float64, orig_data[3:end, n])
end
first(orig_data, 6)

## ------------------------------------------------------------------
# Convert data

# Converting data c(g/L) / MM (g/mol) -> c(mM)
conv_data = DataFrame(orig_data)
conv_data[3:end, :cGLC] .= orig_data[3:end, :cGLC] .* 1e3/ 180;
conv_data[3:end, :sGLC] .= orig_data[3:end, :sGLC] .* 1e3/ 180;
conv_data[3:end, :sSA] .= orig_data[3:end, :sSA] .* 1e3/ 118;
conv_data[3:end, :sAcA] .= orig_data[3:end, :sAcA] .* 1e3/ 60.02;
conv_data[3:end, :sFA] .= orig_data[3:end, :sFA] .* 1e3/ 46.03;
conv_data[3:end, :sMA] .= orig_data[3:end, :sMA] .* 1e3/ 134.09;
conv_data[2, [:cGLC, :sGLC,:sSA,:sAcA,:sFA,:sMA]] .= "mM"
first(conv_data, 6)

## ------------------------------------------------------------------
# Saving

CSV.write(Hd.HEERDEN_CONT_CUL_DATA_CONV_FILE, conv_data, delim = "\t");
println(relpath(Hd.HEERDEN_CONT_CUL_DATA_CONV_FILE), " created!!!")


