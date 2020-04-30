#==
Functions from StatsBase extended to VAR structs
==#

"""
    coef(v::VarEstimate)

Return the matrix of coefficients for VAR model `v`.
"""
coef(v::VarEstimate) = v.B

"""
    stderror(v::VarEstimate)

Return the standard errors for the coefficients of VAR model `v`.
"""
stderror(v::VarEstimate) = v.seB

"""
    residuals(v::VarEstimate)

Return the matrix of residuals for VAR model `v`.
"""
residuals(v::VarEstimate) = v.U

"""
    loglikelihood(v::VarEstimate)

Return the log-likelihood for VAR model `v`.
"""
function loglikelihood(v::VarEstimate)
    U = residuals(v)
    ΣU = v.ΣU
    obs, K = size(U)
    sssr = 0.0
    for t = 1:obs
        sssr += dot(ΣU \ U[t, :], U[t, :])
    end
    ll = -(K * obs / 2) * log(2π) - (obs / 2) * logdet(ΣU) - 0.5 * sssr
end

"""
    aic(v::VarEstimate)

Return Akaike's Information Criterion for VAR model `v`.
"""
function aic(v::VarEstimate)
    obs, K = size(residuals(v))
    nparam = size(coef(v), 1)
    ΣU = v.ΣU
    ΣUml = ((obs - nparam) / obs) * ΣU
    return logdet(ΣUml) + (2 / obs) * nparam * K
end

function coeftable(v::VarEstimate)
    B = coef(v)
    seB = stderror(v)
    m, K = size(B)
    lags = v.lags
    varnames = v.varnames
    obs = v.obs
    dofr = obs - m

    colnms = [
        "Estimate",
        "Std. Error",
        "t value",
        "Pr > |t|",
    ]

    rownms = Vector{String}()
    for l = 1:lags
        for k = 1:K
            push!(rownms, "$(varnames[k]).l$l")
        end
    end
    rownms = v.trend ? ["trend"; rownms] : rownms
    rownms = v.constant ? ["Intercept"; rownms] : rownms

    ctable = Array{CoefTable}(undef, K)
    for k = 1:K
        Bk   = round.(B[:, k], sigdigits=4)
        seBk = round.(seB[:, k], sigdigits=4)
        tk = Bk ./ seBk
        pk = 2 * ccdf.(TDist(dofr), abs.(tk))
        # mat = hcat(Bk, seBk, tk)
        mat = hcat(Bk, seBk, tk, pk)
        # ctable[k] = CoefTable(mat, colnms, rownms)
        ctable[k] = CoefTable(mat, colnms, rownms, 4, 3)
    end
    return ctable
end
