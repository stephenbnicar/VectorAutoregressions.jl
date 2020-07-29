using LinearAlgebra
using StatsBase
using Distributions
cd(@__DIR__)
include("../examples/lutkepohl_example.jl")

h = 4

U = residuals(v)
T, K = size(U)
p = v.lags
nparam = size(v.B, 1)
Z = v.Z
ΣUml = ((T - nparam) / T) * v.ΣU

U2 = vcat(zeros(p + (h - p), K), U)
Ulag = VectorAutoregressions.lag_matrix(U2, h)
ZUlag = [Z Ulag]
AD = (ZUlag' * ZUlag) \ (ZUlag' * U)

E = U - ZUlag * AD
ΣEml = (E' * E) / T

Qlm = T * (K - tr(ΣUml \ ΣEml))
df = h * K^2
pval = ccdf(Chisq(df), Qlm)
