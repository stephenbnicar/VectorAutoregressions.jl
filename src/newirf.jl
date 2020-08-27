# Activate VectorAutoregressions project first
cd(@__DIR__)
include("../examples/lutkepohl_example.jl")
using LinearAlgebra

v = VAR(Y, 2)
p = v.lags
B = v.B
sigma = v.ΣU
K = size(B, 2)
h = 12

phi = VectorAutoregressions.simple_irf(B, K, p, h)
theta = VectorAutoregressions.orthogonal_irf(phi, sigma, K, h)
