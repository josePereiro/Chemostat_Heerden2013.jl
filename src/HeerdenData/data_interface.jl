# This interface tha data stored in data/processed/heerden2013___data/heerden2013___cont_cult_data.tsv
# taken from Heerden2013 https://doi.org/10.1186/1475-2859-12-80.

## ------------------------------------------------------------------
cont_cul_data = DataFrame()
msd_mets = ["GLC", "SA", "AcA", "FA", "MA"]

# Interface
data_ids()::Vector{String} = string.(names(cont_cul_data))
val(id)::Vector{Float64} = parse.(Float64, string.(cont_cul_data[3:end, Symbol(id)]))
val(id, indx::Integer)::Float64 = val(id)[indx]
val(id, indx::Integer, defval)::Float64 = try; return val(id, indx) catch; return defval end
val(id, indx)::Vector{Float64} = val(id)[indx]
name(id)::String = string(cont_cul_data[1, Symbol(id)])
unit(id)::String = string(cont_cul_data[2, Symbol(id)])

# Load data
function _load_cont_cul_data()

    !isfile(HEERDEN_CONT_CUL_DATA_CONV_FILE) && return cont_cul_data
    
    ## ------------------------------------------------------------------
    global cont_cul_data = CSV.read(HEERDEN_CONT_CUL_DATA_CONV_FILE, DataFrame; delim = '\t', comment = "#")

    ## ------------------------------------------------------------------
    # include xi
    cont_cul_data[!, :xi] = ["Cell specific dilution rate"; "gDW/ L hr"; val(:DCW) ./ val(:D)];

    # Deleting data set with NaN xi values
    cont_cul_data = cont_cul_data[map(x -> x isa AbstractString || !isnan(x), cont_cul_data.xi),:]

    ## ------------------------------------------------------------------
    # Converting to flux u = (s - c)/ξ (+) means uptake
    let
        xi = val(:xi)
        for met in msd_mets
            umet = Symbol("u$met")
            cmet = Symbol("c$met")
            smet = Symbol("s$met")
            s = val(smet)
            c = string(cmet) in names(cont_cul_data) ? val(cmet) : zero(s)
            cont_cul_data[umet] = ["Exchange"; "mmol/gDW hr"; (s .- c) ./ xi];
        end
    end

    return cont_cul_data
end
