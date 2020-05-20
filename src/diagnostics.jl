"""
    sic(v::VarEstimate)

Return the Schwarz (Bayesian) Information Criterion for VAR model `v`.
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

Return the Hannan-Quinn Criterion for VAR model `v`.
"""
function hqc(v::VarEstimate)
    obs, K = size(residuals(v))
    nparam = size(coef(v), 1)
    ΣU = v.ΣU
    ΣUml = ((obs - nparam) / obs) * ΣU
    return logdet(ΣUml) + (2 * log(log(obs)) / obs) * nparam * K
end

"""
    VarStabilityCheck

# Fields
- `isstable::Bool`
- `eigenvals::Array{Number}`
- `eigenmod::Array{Float64}`
"""
struct VarStabilityCheck
    isstable::Bool
    eigenvals::Array{Number}
    eigenmod::Array{Float64}
end

"""
    checkstable(v::VarEstimate) -> VarStabilityCheck

Check the stability of the estimated VAR model `v`.
"""
function checkstable(v::VarEstimate)
    # See Lutkepohl (2006) p.15ff
    p = v.lags
    B = coef(v)
    nparam, K = size(B)
    startidx = nparam - K * p + 1
    A = B[startidx:end, :]
    Acomp = companion(A, p)
    eigenvals = eigvals(Acomp)
    eigenmod = abs.(eigenvals)
    isstable = all(eigenmod .< 1)
    VarStabilityCheck(isstable, eigenvals, eigenmod)
end

function show(io::IO, obj::VarStabilityCheck)
    # K = length(stab.ynames)
    E = obj.eigenvals
    Emod = round.(unique(obj.eigenmod), digits = 3)
    println(io, typeof(obj))
    if obj.isstable
        println(io, "VAR is stable")
    else
        println(io, "VAR is not stable")
    end
    println(io, "Modulus of eigenvalues (dropping repeated values):")
    println(io, sort(Emod))
end

"""
    portmanteau_test(v::VarEstimate, h::Int)

Conduct a multivariate portmanteau test for residual autocorrelation up to lag
    `h` for VAR model `v`.
Implements the adjusted test described on p.171 of Lutkepohl (2006).
"""
function portmanteau_test(v::VarEstimate, h::Int)
    U = residuals(v)
    T, K = size(U)
    p = v.lags
    C = permutedims(crosscov(U, U, collect(0:h), demean = false), [3, 2, 1])
    Q = 0.0
    for i = 1:h
        Q += tr(C[:, :, i+1]' * (C[:, :, 1] \ C[:, :, i+1] / C[:, :, 1])) * (T^2 / (T - i))
    end
    df = K^2 * (h - p)
    pval = ccdf(Chisq(df), Q)
    return Q, df, pval
end
