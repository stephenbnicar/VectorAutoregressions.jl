function loglikelihood(v::VAR)
    U = v.residuals
    Σᵤ = v.vcov_residuals
    obs, K = size(U)
    sssr = 0.0
    for t = 1:obs
        sssr += dot(Σᵤ \ U[t, :], U[t, :])
    end
    ll = -(K * obs / 2) * log(2π) - (obs / 2) * logdet(Σᵤ) - 0.5 * sssr
end

function aic(v::VAR)
    obs, K = size(v.residuals)
    nparam = size(v.coef, 1)
    Σᵤ = v.vcov_residuals
    Σ = ((obs - nparam) / obs) * Σᵤ
    return logdet(Σ) + (2 / obs) * nparam * K
end

function bic(v::VAR)
    obs, K = size(v.residuals)
    nparam = size(v.coef, 1)
    Σᵤ = v.vcov_residuals
    Σ = ((obs - nparam) / obs) * Σᵤ
    return logdet(Σ) + (log(obs) / obs) * nparam * K
end

function hqc(v::VAR)
    obs, K = size(v.residuals)
    nparam = size(v.coef, 1)
    Σᵤ = v.vcov_residuals
    Σ = ((obs - nparam) / obs) * Σᵤ
    return logdet(Σ) + (2 * log(log(obs)) / obs) * nparam * K
end
