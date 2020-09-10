"""
    sic(v::VAREstimate)

Return the Schwarz (Bayesian) Information Criterion for VAR model `v`.
"""
function sic(v::VAREstimate)
    obs, K = size(residuals(v))
    nparam = size(coef(v), 1)
    ΣU = v.ΣU
    ΣUml = ((obs - nparam) / obs) * ΣU
    return logdet(ΣUml) + (log(obs) / obs) * nparam * K
end

"""
    hqc(v::VAREstimate)

Return the Hannan-Quinn Criterion for VAR model `v`.
"""
function hqc(v::VAREstimate)
    obs, K = size(residuals(v))
    nparam = size(coef(v), 1)
    ΣU = v.ΣU
    ΣUml = ((obs - nparam) / obs) * ΣU
    return logdet(ΣUml) + (2 * log(log(obs)) / obs) * nparam * K
end

"""
    StabilityCheck(v::VAREstimate) -> StabilityCheck

Check the stability of the estimated VAR model `v`.

# Fields
- `isstable::Bool`
- `eigenvals::Array{Number}`
- `eigenmod::Array{Float64}`
"""
struct StabilityCheck
    isstable::Bool
    eigenvals::Array{Number}
    eigenmod::Array{Float64}
end

function StabilityCheck(v::VAREstimate)
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
    StabilityCheck(isstable, eigenvals, eigenmod)
end

function show(io::IO, obj::StabilityCheck)
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
    PortmanteauTest(v::VAREstimate, h::Int) -> PortmanteauTest

Conduct a multivariate portmanteau test for residual autocorrelation up to lag
    `h` for VAR model `v`.
Implements the adjusted test described on p.171 of Lütkepohl (2006).

# Fields
- `U::Matrix`
- `h::Int`
- `Q::Float64`
- `df::Int`
- `pval::Float64`
"""
struct PortmanteauTest
    U::Matrix
    h::Int
    Q::Float64
    df::Int
    pval::Float64
end

function PortmanteauTest(v::VAREstimate, h::Int)
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
    PortmanteauTest(U, h, Q, df, pval)
end

function show(io::IO, obj::PortmanteauTest)
    println(io, typeof(obj))
    println(io, "Multivariate test for residual autocorrelation.")
    println(io, "- Null is no autocorrelation")
    println(io, "- Lags: $(obj.h)")
    println(io, "- Q stat: $(round(obj.Q; digits=3))")
    println(io, "- p-value: $(round(obj.pval; digits=3))")
    outcome = obj.pval > 0.05 ? "Fail to reject" : "Reject"
    println(io, "$outcome the null at the 5% level")
end

"""
    LMCorrTest(v::VAREstimate, h::Int; smallsample = false) -> LMCorrTest

Conduct a multivariate LM test (Breusch-Godfrey test) for residual autocorrelation
    up to lag `h` for VAR model `v`.
If `smallsample = true`, the alternative test statistic detailed in Lütkepohl (2006),
    p.173 is used.

# Fields
- `h::Int`
- `Q::Float64`
- `df1::Any`
- `df2::Any`
- `pval::Float64`
- `smallsample::Bool`
"""
struct LMCorrTest
    h::Int
    Q::Float64
    df1::Any
    df2::Any
    pval::Float64
    smallsample::Bool
end

function LMCorrTest(v::VAREstimate, h::Int; smallsample = false)
    U = residuals(v)
    T, K = size(U)
    p = v.lags
    nparam = size(v.B, 1)
    Z = v.Z

    U2 = vcat(zeros(p + (h - p), K), U)
    Ulag = VectorAutoregressions.lag_matrix(U2, h)
    ZUlag = [Z Ulag]
    AD = (ZUlag' * ZUlag) \ (ZUlag' * U)
    E = U - ZUlag * AD

    ΣUml = ((T - nparam) / T) * v.ΣU
    ΣEml = (E' * E) / T

    if smallsample
        # Test stat from Lütkepohl 2006, p.173
        s = sqrt((K^4 * h^2 - 4) / (K^2 + (K^2 * h^2) - 5))
        N = T - (K*p) - 1 - (K*h) - (K - K*h +1)/2
        df1 = h * K^2
        df2 = (N * s) - ((K^2 * h) / 2) + 1
        Q = ((det(ΣUml) / det(ΣEml))^(1 / s) - 1) * (df2 / df1)
        pval = ccdf(FDist(df1, df2), Q)
    else
        # Test stat from Killian and Lütkepohl 2017, p.54
        Q = T * (K - tr(ΣUml \ ΣEml))
        df1 = h * K^2
        df2 = nothing
        pval = ccdf(Chisq(df1), Q)
    end
    return LMCorrTest(h, Q, df1, df2, pval, smallsample)
end

function show(io::IO, obj::LMCorrTest)
    println(io, typeof(obj))
    println(io, "Multivariate LM test for residual autocorrelation.")
    println(io, "- Null is no autocorrelation")
    println(io, "- Lags: $(obj.h)")
    println(io, "- Q stat: $(round(obj.Q; digits=3))")
    println(io, "- p-value: $(round(obj.pval; digits=3))")
    if obj.smallsample
        println(io, "- Small sample test")
    end
    outcome = obj.pval > 0.05 ? "Fail to reject" : "Reject"
    println(io, "$outcome the null at the 5% level")
end
