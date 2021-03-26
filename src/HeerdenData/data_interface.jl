# This interface tha data stored in data/processed/heerden2013___data/heerden2013___cont_cult_data.tsv
# taken from Heerden2013 https://doi.org/10.1186/1475-2859-12-80.

## ------------------------------------------------------------------
CONT_CUL_DATA = DataFrame()
MSD_METS = ["GLC", "SUCC", "AC", "FORM", "MAL"]
EXPS = 1:13

# Interface
data_ids()::Vector{String} = string.(names(CONT_CUL_DATA))
val(id)::Vector{Float64} = parse.(Float64, string.(CONT_CUL_DATA[3:end, Symbol(id)]))
val(id, indx::Integer)::Float64 = val(id)[indx]
val(id, indx::Integer, defval)::Float64 = try; return val(id, indx) catch; return defval end
val(id, indx)::Vector{Float64} = val(id)[indx]
name(id)::String = string(CONT_CUL_DATA[1, Symbol(id)])
unit(id)::String = string(CONT_CUL_DATA[2, Symbol(id)])


for preffix in ["u", "s", "c"]
    fname = Symbol(preffix, :val)
    @eval $fname(id, args...) = val(string($preffix, id), args...)
end

# Load data
function _load_cont_cul_data()

    !isfile(HEERDEN_CONT_CUL_DATA_CONV_FILE) && return CONT_CUL_DATA
    
    ## ------------------------------------------------------------------
    global CONT_CUL_DATA = CSV.read(HEERDEN_CONT_CUL_DATA_CONV_FILE, DataFrame; delim = '\t', comment = "#")

    ## ------------------------------------------------------------------
    # include xi
    CONT_CUL_DATA[!, :xi] = ["Cell specific dilution rate"; "gDW/ L hr"; val(:DCW) ./ val(:D)];

    # Deleting data set with NaN xi values
    CONT_CUL_DATA = CONT_CUL_DATA[map(x -> x isa AbstractString || !isnan(x), CONT_CUL_DATA.xi),:]

    ## ------------------------------------------------------------------
    # Converting to flux u = (s - c)/ξ (+) means uptake
    let
        xi = val(:xi)
        for met in MSD_METS
            umet = Symbol("u$met")
            smet = Symbol("s$met")
            cmet = Symbol("c$met")
            s = val(smet)
            c = string(cmet) in names(CONT_CUL_DATA) ? val(cmet) : zero(s)
            CONT_CUL_DATA[!, umet] = ["Exchange"; "mmol/gDW hr"; (s .- c) ./ xi];
        end
    end

    return CONT_CUL_DATA
end