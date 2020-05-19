"""
    LagSelectionCriteria

# Fields
- `maxlag::Int`
- `table::DataFrame`
- `selection::Dict`
"""
struct LagSelectionCriteria
    maxlag::Int
    table::DataFrame
    selection::Dict
end

"""
    lagselect(data, maxlag; constant = true, trend = false) -> LagSelectionCriteria

Calculate AIC, SIC, and HQC lag selection criteria for an unrestricted VAR.

# Arguments
- `data` : `DataFrame` or `TimeArray` of observations on endogenous variables
- `maxlag::Int` : the maximum number of lags to estimate
- `constant::Bool = true` : include an intercept term
- `trend::Bool = false` : include a linear trend
"""
function lagselect(data::DataFrame, maxlag::Int; constant::Bool = true, trend::Bool = false)
    if maxlag < 1
        error("maxlag must be ≥ 1")
    end
    AIC = zeros(maxlag)
    HQC = zeros(maxlag)
    SIC = zeros(maxlag)
    for lag = 1:maxlag
        dmat = data[(maxlag-lag+1):end, :]
        v = VAR(dmat, lag; constant = constant, trend = trend)
        AIC[lag] = aic(v)
        HQC[lag] = hqc(v)
        SIC[lag] = sic(v)
    end
    table = DataFrame(lag = collect(1:maxlag), AIC = AIC, HQC = HQC, SIC = SIC)
    selection = Dict(cn => argmin(table[!, cn]) for cn in names(table)[2:end])
    LagSelectionCriteria(maxlag, table, selection)
end

function lagselect(data::TimeArray, maxlag::Int; constant::Bool = true, trend::Bool = false)
    data = DataFrame(data)[:, 2:end]
    lagselect(data, maxlag; constant = constant, trend = trend)
end

function show(io::IO, ls::LagSelectionCriteria)
    println(io, "VAR Lag Selection")
    println(io, "-------------------------")
    println(io, "maxlag: ", ls.maxlag)
    println(io, "Criterion   Lag Selection")
    println(io, "-------------------------")
    for k in keys(ls.selection)
        println(io, rpad(string(k), 16), ls.selection[k])
    end
    println(io, "-------------------------")
end
