# function criteria()
#     tobs, K = size(v.data)
#     p = v.lags
#     nobs = tobs - p
#     nparam = (K*p) + v.constant + v.trend
#     Σ_U = v.vcov_residuals
#     # MLE estimate of SigmaU
#     Σ = ((nobs - nparam)/nobs)*Σ_U
#
#     FPE  = det(Σ)*(((nobs + nparam)/(nobs - nparam))^K)*1e11
#     AIC  = logdet(Σ) + (2/nobs)*nparam*K
#     HQIC = logdet(Σ) + (2log(log(nobs))/nobs)*nparam*K
#     SIC  = logdet(Σ) + (log(nobs)/nobs)*nparam*K
#
#     return Dict("FPE" => FPE, "AIC" => AIC, "SIC" => SIC, "HQIC" => HQIC)
# end
