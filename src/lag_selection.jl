"""
    LagSelection(data, maxlag; constant = true, trend = false) -> LagSelection

Calculate AIC, SIC, and HQC lag selection criteria for an unrestricted VAR.

# Arguments
- `data` : `DataFrame` or `TimeArray` of observations on endogenous variables
- `maxlag::Int` : the maximum number of lags to estimate
- `constant::Bool = true` : include an intercept term
- `trend::Bool = false` : include a linear trend

# Fields
- `maxlag::Int`
- `table::DataFrame`
- `selection::Dict`
"""
struct LagSelection
    maxlag::Int
    table::DataFrame
    selection::Dict
end

function LagSelection(data::DataFrame, maxlag::Int; constant::Bool = true, trend::Bool = false)
    if maxlag < 1
        error("maxlag must be ≥ 1")
    end

    table, selection = lagselect(data, maxlag, constant, trend)
    LagSelection(maxlag, table, selection)
end

function LagSelection(data::TimeArray, maxlag::Int; constant::Bool = true, trend::Bool = false)
    data = DataFrame(data)[:, 2:end]
    LagSelection(data, maxlag; constant = constant, trend = trend)
end

function lagselect(data, maxlag, constant, trend)
    critname = ["AIC", "HQIC", "SIC"]
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
    return table, selection
end

function show(io::IO, ls::LagSelection)
    println(io, typeof(ls))
    println(io, "Maximum lags: ", ls.maxlag)
    for k in keys(ls.selection)
        print(io, k, ": ", rpad(ls.selection[k], 3))
    end
    print(io, '\n')
end
