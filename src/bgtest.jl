using LinearAlgebra
using StatsBase
cd(@__DIR__)
include("../examples/lutkepohl_example.jl")

U = residuals(v)
T = v.obs
h = 2

function makeF(T, h)
    F = Matrix{Float64}(I, T, T)
    for i in 1:h
        Fi = [zeros(i, T); Matrix{Float64}(I, (T-i), (T-i)) zeros(T-i, i)]
        F = hcat(F, Fi)
    end
    F[:, (T+1):end]
end

F = makeF(T, h)

# See (4.4.2) on p.158; that formula is missing 1/T # hide
C = U' * F * kron(Matrix{Float64}(I, h, h), U) / T

Calt = permutedims(crosscov(U, U, collect(0:h), demean = false), [3, 2, 1])
