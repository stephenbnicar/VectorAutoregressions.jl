#==
Functions from StatsBase extended to VAR structs
==#

"""
    coef(v::VAREstimate)

Return the matrix of coefficients for VAR model `v`.
"""
coef(v::VAREstimate) = v.B

"""
    stderror(v::VAREstimate)

Return the standard errors for the coefficients of VAR model `v`.
"""
stderror(v::VAREstimate) = v.seB

"""
    residuals(v::VAREstimate)

Return the matrix of residuals for VAR model `v`.
"""
residuals(v::VAREstimate) = v.U

"""
    fitted(v::VAREstimate)

Return the fitted values for VAR model `v`.
"""
fitted(v::VAREstimate) = v.Yhat

"""
    loglikelihood(v::VAREstimate)

Return the log-likelihood for VAR model `v`.
"""
function loglikelihood(v::VAREstimate)
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
    aic(v::VAREstimate)

Return Akaike's Information Criterion for VAR model `v`.
"""
function aic(v::VAREstimate)
    obs, K = size(residuals(v))
    nparam = size(coef(v), 1)
    ΣU = v.ΣU
    ΣUml = ((obs - nparam) / obs) * ΣU
    return logdet(ΣUml) + (2 / obs) * nparam * K
end

function coeftable(v::VAREstimate)
    B = coef(v)
    seB = stderror(v)
    m, K = size(B)
    lags = v.lags
    ynames = v.ynames
    xnames = v.xnames
    obs = v.obs
    dofr = obs - m

    colnms = ["Estimate", "Std. Error", "t value", "Pr > |t|"]

    rownms = Vector{String}()
    for l = 1:lags
        for k = 1:K
            push!(rownms, "$(ynames[k]).l$l")
        end
    end
    rownms = !isa(v.X, Nothing) ? [xnames; rownms] : rownms
    rownms = v.trend ? ["trend"; rownms] : rownms
    rownms = v.constant ? ["intercept"; rownms] : rownms

    ctable = Array{CoefTable}(undef, K)
    for k = 1:K
        Bk = round.(B[:, k], digits = 4)
        seBk = round.(seB[:, k], digits = 4)
        tk = Bk ./ seBk
        pk = 2 * ccdf.(TDist(dofr), abs.(tk))
        # mat = hcat(Bk, seBk, tk)
        mat = hcat(Bk, seBk, tk, pk)
        # ctable[k] = CoefTable(mat, colnms, rownms)
        ctable[k] = CoefTable(mat, colnms, rownms, 4, 3)
    end
    return ctable
end
