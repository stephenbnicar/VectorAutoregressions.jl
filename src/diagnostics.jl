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

`struct` to hold the results of [`checkstable`](@ref).

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
    PortmanteauTest

`struct` to hold the results of [`portmanteau_test`](@ref).

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

"""
    portmanteau_test(v::VarEstimate, h::Int) -> PortmanteauTest

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
    return PortmanteauTest(U, h, Q, df, pval)
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

struct BreuschGodfreyTest
    U::Matrix
    h::Int
    Q::Float64
    df::Int
    pval::Float64
end

function bg_test(v::VarEstimate, h::Int)
    U = residuals(v)
    T, K = size(U)
    p = v.lags
    nparam = size(v.B, 1)
    Z = v.Z

    U2 = vcat(zeros(p+(h-p), K), U)
    Ulag = VectorAutoregressions.lag_matrix(U2, h)
    ZUlag = [Z Ulag]
    AD = (ZUlag' * ZUlag) \ (ZUlag' * U)
    E = U - ZUlag * AD

    ΣUml = ((T - nparam) / T) * v.ΣU
    ΣEml = (E' * E) / T

    Q = T * (K - tr(ΣUml \ ΣEml))
    df = h * K^2
    pval = ccdf(Chisq(df), Q)
    return BreuschGodfreyTest(U, h, Q, df, pval)
end

function show(io::IO, obj::BreuschGodfreyTest)
    println(io, typeof(obj))
    println(io, "Multivariate test for residual autocorrelation.")
    println(io, "- Null is no autocorrelation")
    println(io, "- Lags: $(obj.h)")
    println(io, "- Q stat: $(round(obj.Q; digits=3))")
    println(io, "- p-value: $(round(obj.pval; digits=3))")
    outcome = obj.pval > 0.05 ? "Fail to reject" : "Reject"
    println(io, "$outcome the null at the 5% level")
end
