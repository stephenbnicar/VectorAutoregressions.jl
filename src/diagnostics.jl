"""
    loglikelihood(v::VectorAutoregressions.VAR)

Calculate log likelihood for the VAR model `v`.
"""
function loglikelihood(v::VAR)
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
    aic(v::VectorAutoregressions.VAR)

Calculate Akaike's Information Criterion for the VAR model `v`.
"""
function aic(v::VAR)
    obs, K = size(residuals(v))
    nparam = size(coef(v), 1)
    ΣU = v.ΣU
    ΣUml = ((obs - nparam) / obs) * ΣU
    return logdet(ΣUml) + (2 / obs) * nparam * K
end

"""
    sic(v::VectorAutoregressions.VAR)

Calculate the Schwarz (Bayesian) Information Criterion for the VAR model `v`.
"""
function sic(v::VAR)
    obs, K = size(residuals(v))
    nparam = size(coef(v), 1)
    ΣU = v.ΣU
    ΣUml = ((obs - nparam) / obs) * ΣU
    return logdet(ΣUml) + (log(obs) / obs) * nparam * K
end

"""
    hqc(v::VectorAutoregressions.VAR)

Calculate the Hannan-Quinn Criterion for the VAR model `v`.
"""
function hqc(v::VAR)
    obs, K = size(residuals(v))
    nparam = size(coef(v), 1)
    ΣU = v.ΣU
    ΣUml = ((obs - nparam) / obs) * ΣU
    return logdet(ΣUml) + (2 * log(log(obs)) / obs) * nparam * K
end
