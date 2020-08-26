# Activate VectorAutoregressions project first
cd(@__DIR__)
include("../examples/lutkepohl_example.jl")
using LinearAlgebra

v = VAR(Y, 2)
p = v.lags
B = v.B
K = size(B, 2)
h = 3
