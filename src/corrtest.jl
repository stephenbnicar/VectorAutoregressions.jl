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
