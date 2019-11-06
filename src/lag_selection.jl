struct LagSelectionCriteria
    table::DataFrame
    selection::Dict
end

function lagselect(Y, maxlag; constant = true, trend = false)
    critname = ["AIC", "HQC", "SIC"]
    AIC = []
    HQC = []
    SIC = []
    for lag = 1:maxlag
        dmat = Y[(maxlag+1-lag):end, :]
        v = VAR(dmat, lag; constant = constant, trend = trend)
        infcrit = criteria(v)
        push!(AIC, infcrit["AIC"])
        push!(HQC, infcrit["HQC"])
        push!(SIC, infcrit["SIC"])
    end
    table = DataFrame(lag = collect(1:maxlag), AIC = AIC, HQC = HQC, SIC = SIC)
    selection = Dict(cn => argmin(table[cn]) for cn in names(table)[2:end])
    return LagSelectionCriteria(table, selection)
end

function criteria(v::VAR)
    tobs, K = size(v.data)
    p = v.lags
    obs = tobs - p
    nparam = (K * p) + v.constant + v.trend
    Σᵤ = v.vcov_resid
    # MLE estimate of Σᵤ
    Σ = ((obs - nparam) / obs) * Σᵤ

    AIC = logdet(Σ) + (2 / obs) * nparam * K
    HQC = logdet(Σ) + (2 * log(log(obs)) / obs) * nparam * K
    SIC = logdet(Σ) + (log(obs) / obs) * nparam * K

    return Dict("AIC" => AIC, "SIC" => SIC, "HQC" => HQC)
end

function show(io::IO, ls::LagSelectionCriteria)
    println(io, "VAR Lag Selection Criteria")
    println(io, "---------------------------")
    println(io, "Criteria   Lag")
    println(io, "---------------------------")
    for k in keys(ls.selection)
        println(io, rpad(string(k), 11), ls.selection[k])
    end
    println(io, "---------------------------")
end
