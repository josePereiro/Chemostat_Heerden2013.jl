# This interface tha data stored in data/processed/heerden2013___data/heerden2013___cont_cult_data.tsv
# taken from Heerden2013 https://doi.org/10.1186/1475-2859-12-80.

# ---
# ## Load data
# ---
cont_cul_data = CSV.read(HEERDEN_CONT_CUL_DATA_CONV_FILE, DataFrame; delim = '\t', comment = "#")

# ---
# ## Interface
# ---
msd_mets = ["GLC", "SA", "AcA", "FA", "MA"]
data_ids()::Vector{String} = string.(names(cont_cul_data))
val(id)::Vector{Float64} = parse.(Float64, string.(cont_cul_data[3:end, Symbol(id)]))
val(id, indx::Integer)::Float64 = val(id)[indx]
val(id, indx)::Vector{Float64} = val(id)[indx]
name(id)::String = string(cont_cul_data[1, Symbol(id)])
unit(id)::String = string(cont_cul_data[2, Symbol(id)])

# include xi
cont_cul_data[!, :xi] = ["Cell specific dilution rate"; "gDW/ L hr"; val(:DCW) ./ val(:D)];

# Deleting data set with NaN xi values
cont_cul_data = cont_cul_data[map(x -> x isa AbstractString || !isnan(x), cont_cul_data.xi),:]