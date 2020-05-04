cd(@__DIR__)
include("../examples/kilian_example.jl")

using LinearAlgebra
using StatsBase

U = residuals(v)
T, K = size(U)

# See Lutkepohl p.158, 161
C1 = Matrix{Float64}(undef, K, K)
C1[1,1] = dot(U[1:T-1,1], U[2:T,1])/T
C1[2,1] = dot(U[1:T-1,1], U[2:T,2])/T
C1[3,1] = dot(U[1:T-1,1], U[2:T,3])/T

X = dropdims(permutedims(crosscov(U, U, [1], demean = false), [3, 2, 1]), dims=3)

function portmanteau_test1(v::VarEstimate, h)
    U = residuals(v)
    T, K = size(U)
    C = Array{Float64}(undef, K, K, h+1)
    C[:, :, 1] = (U'*U)/T
    # invC0 = inv(C0)
    # Q = 0.0
    for i in 1:h
        F = [zeros(i, T); Matrix(1.0I, (T-i), (T-i)) zeros(T-i, i)]
        C[:, :, i+1] = (U'*F*U)/T
        # Q += tr(Ci'*invC0*Ci*invC0)*(T^2/(T - i))
    end
    # df = K^2*(h-p)
    # pval = 1 - cdf(Chisq(df), Q)
    # return Q, df, pval
    return C
end

# This is faster, fewer allocations
function portmanteau_test2(v::VarEstimate, h)
    U = residuals(v)
    T, K = size(U)
    C = permutedims(crosscov(U, U, collect(0:h), demean = false), [3, 2, 1])
end
