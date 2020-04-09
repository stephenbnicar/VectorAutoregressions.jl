"""
    sic(v::VarEstimate)

Calculate the Schwarz (Bayesian) Information Criterion for the VAR model `v`.
"""
function sic(v::VarEstimate)
    obs, K = size(residuals(v))
    nparam = size(coef(v), 1)
    ΣU = v.ΣU
    ΣUml = ((obs - nparam) / obs) * ΣU
    return logdet(ΣUml) + (log(obs) / obs) * nparam * K
end

"""
    hqc(v::VarEstimate)

Calculate the Hannan-Quinn Criterion for the VAR model `v`.
"""
function hqc(v::VarEstimate)
    obs, K = size(residuals(v))
    nparam = size(coef(v), 1)
    ΣU = v.ΣU
    ΣUml = ((obs - nparam) / obs) * ΣU
    return logdet(ΣUml) + (2 * log(log(obs)) / obs) * nparam * K
end
