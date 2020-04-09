coef(v::VarEstimate) = v.B
stderror(v::VarEstimate) = v.seB
residuals(v::VarEstimate) = v.U

"""
    loglikelihood(v::VarEstimate)

Calculate log likelihood for the VAR model `v`.
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

Calculate Akaike's Information Criterion for the VAR model `v`.
"""
function aic(v::VarEstimate)
    obs, K = size(residuals(v))
    nparam = size(coef(v), 1)
    ΣU = v.ΣU
    ΣUml = ((obs - nparam) / obs) * ΣU
    return logdet(ΣUml) + (2 / obs) * nparam * K
end
