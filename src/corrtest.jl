cd(@__DIR__)
include("../examples/lutkepohl_example.jl")

using LinearAlgebra
using StatsBase
using Distributions

U = residuals(v)
T, K = size(U)

# See Lutkepohl p.158, 161
function portmanteau_test(v::VarEstimate, h)
    U = residuals(v)
    T, K = size(U)
    p = v.lags
    C = permutedims(crosscov(U, U, collect(0:h), demean = false), [3, 2, 1])
    Q = 0.0
    for i in 1:h
        Q += tr(C[:,:,i+1]'*(C[:,:,1]\C[:,:,i+1]/C[:,:,1]))*(T^2/(T - i))
    end
    df = K^2*(h-p)
    pval = ccdf(Chisq(df), Q)
    return Q, df, pval
end
