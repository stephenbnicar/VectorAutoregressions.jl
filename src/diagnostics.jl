function loglikelihood(residuals, vcov_resid)
    # U = v.residuals
    # Σ = v.vcov_residuals
    uobs, K = size(residuals)
    sssr = 0.0
    for t in 1:uobs
        sssr += dot(vcov_resid\residuals[t, :], residuals[t, :])
    end
    ll = -(K*uobs/2)*log(2π) - (uobs/2)*logdet(vcov_resid) - 0.5*sssr
end

# function criteria()
#     tobs, K = size(v.data)
#     p = v.lags
#     nobs = tobs - p
#     nparam = (K*p) + v.constant + v.trend
#     Σ_U = v.vcov_residuals
#     # MLE estimate of SigmaU
#     Σ = ((nobs - nparam)/nobs)*Σ_U
#
#     AIC  = logdet(Σ) + (2/nobs)*nparam*K
#     HQIC = logdet(Σ) + (2log(log(nobs))/nobs)*nparam*K
#     SIC  = logdet(Σ) + (log(nobs)/nobs)*nparam*K
#
#     return Dict("FPE" => FPE, "AIC" => AIC, "SIC" => SIC, "HQIC" => HQIC)
# end
