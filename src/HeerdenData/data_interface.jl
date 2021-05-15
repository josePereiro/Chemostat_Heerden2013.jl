# This interface tha data stored in data/processed/heerden2013___data/heerden2013___cont_cult_data.tsv
# taken from Heerden2013 https://doi.org/10.1186/1475-2859-12-80.

## ------------------------------------------------------------------
CONT_CUL_DATA = DataFrame()
# msd_mets = ["GLC", "SUCC", "AC", "FORM", "MAL"]
msd_mets = ["GLC", "SUCC", "AC", "FORM"]
EXPS = 1:13

# Interface
parse_id(id) = string(id) == "X" ? :DCW : Symbol(id)
data_ids()::Vector{String} = string.(names(CONT_CUL_DATA))
val(id)::Vector{Float64} = parse.(Float64, string.(CONT_CUL_DATA[3:end, parse_id(id)]))
val(id, indx::Integer)::Float64 = val(id)[indx]
val(id, indx::Integer, defval)::Float64 = try; return val(id, indx) catch; return defval end
val(id, indx)::Vector{Float64} = val(id)[indx]
name(id)::String = string(CONT_CUL_DATA[1, parse_id(id)])
unit(id)::String = string(CONT_CUL_DATA[2, parse_id(id)])

for preffix in ["u", "s", "c"]
    fname = Symbol(preffix, :val)
    @eval $fname(id, args...) = val(string($preffix, id), args...)
end

ciD_X(id) = [cval(id, exp, 0.0) * val(:D, exp) / val(:DCW, exp) for exp in EXPS]
ciD_X(id, exp) = cval(id, exp, 0.0) * val(:D, exp) / val(:DCW, exp)

# Load data
function _load_cont_cul_data()
    
    CONT_CUL_DATA_CONV_FILE = procdir("heerden2013___cont_cult_data.tsv")
    !isfile(CONT_CUL_DATA_CONV_FILE) && return CONT_CUL_DATA
    
    ## ------------------------------------------------------------------
    global CONT_CUL_DATA = CSV.read(CONT_CUL_DATA_CONV_FILE, DataFrame; 
        delim = '\t', comment = "#"
    )
    
    ## ------------------------------------------------------------------
    # include xi
    CONT_CUL_DATA[!, :xi] = [
        "Cell specific dilution rate"; 
        "gDW/ L hr"; 
        val(:DCW) ./ val(:D)
    ];
    
    ## ------------------------------------------------------------------
    # Deleting data set with NaN xi values
    CONT_CUL_DATA = CONT_CUL_DATA[map(x -> x isa AbstractString || !isnan(x), CONT_CUL_DATA.xi),:]

    ## ------------------------------------------------------------------
    # Converting to flux u = (s - c)/ξ (+) means uptake
    let
        xi = val(:xi)
        for met in msd_mets
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