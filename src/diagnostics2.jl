"""
    ResidualCorrelationTest(v::VAR)

Calculate the test statistic and p-value for a multivariate portmanteau test for residual autocorrelation.
Implements the adjusted test described on p.171 of Lutkepohl (2006).
"""
struct ResidualCorrelationTest
    h::Int
    Q::Float64
    pvalQ::Float64
end

function ResidualCorrelationTest(v::T, h::Int) where {T <: AbstractVARModel}
    p = v.lags
    U = v.residuals
    uobs, K = size(U)
    Q, dfQ, pvalQ = portmanteau_test(U, K, p, uobs, h)
    ResidualCorrelationTest(h, Q, pvalQ)
end

function portmanteau_test(U, K, p, uobs, h)
    C0 = (U'*U)/uobs
    invC0 = inv(C0)
    Q = 0.0
    for i in 1:h
        F = [zeros(i, uobs); Matrix{Float64}(I, (uobs-i), (uobs-i)) zeros(uobs-i, i)]
        Ci = (U'*F*U)/uobs
        Q += tr(Ci'*invC0*Ci*invC0)*(uobs^2/(uobs - i))
    end
    df = K^2*(h-p)
    pval = 1 - cdf(Chisq(df), Q)
    return Q, df, pval
end

function show(io::IO, rct::ResidualCorrelationTest)
    println(io, "VAR Residual Correlation Test")
    println(io, "--------------------------------")
    println(io, "Multivariate portmanteau test")
    println(io, "(Null is no autocorrelation)")
    println(io, "--------------------------------")
    println(io, "Q statistic:", lpad(string(round(rct.Q; digits=3)), 12))
    println(io, "P-value:    ", lpad(string(round(rct.pvalQ; digits=3)), 12))
    outcome = rct.pvalQ > 0.05 ? "Fail to reject" : "Reject"
    println(io, "$outcome the null at the 5% level")
    println(io, "--------------------------------")
end


struct ResidualNormalityTest
    λ_s::Float64
    λ_k::Float64
    λ_sk::Float64
    pval1::Float64
    pval2::Float64
    pval3::Float64
end

function ResidualNormalityTest(v::T) where {T <: AbstractVARModel}
    U = v.residuals
    Σ_U = v.vcov_residuals
    λ_s, λ_k, λ_sk, pval1, pval2, pval3 = multivariate_normality_test(U, Σ_U)
    ResidualNormalityTest(λ_s, λ_k, λ_sk, pval1, pval2, pval3)
end

function multivariate_normality_test(U, Σ)
    # see Lutkepohl pp. 177-181
    uobs, K = size(U)
    P   = cholesky(Σ).L
    W   = transpose(P\U')
    W3  = W.^3
    W4  = W.^4
    b1  = dropdims(sum(W3, dims=1)/uobs, dims=1)
    b2  = dropdims(sum(W4, dims=1)/uobs, dims=1)
    λ_s = dot(b1, b1)*(uobs/6)
    λ_k = dot(b2 .- 3, b2 .- 3)*(uobs/24)
    λ_sk = λ_s + λ_k

    pval1 = 1 - cdf(Chisq(K), λ_s)
    pval2 = 1 - cdf(Chisq(K), λ_k)
    pval3 = 1 - cdf(Chisq(2K), λ_sk)

    return λ_s, λ_k, λ_sk, pval1, pval2, pval3
end
