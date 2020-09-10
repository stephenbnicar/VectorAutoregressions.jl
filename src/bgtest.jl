using LinearAlgebra
using StatsBase
using Distributions
cd(@__DIR__)
include("../examples/lutkepohl_example.jl")

h = 8

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

# Test stat from Lütkepohl 2006, p.173
s = sqrt((K^4 * h^2 - 4) / (K^2 + (K^2 * h^2) - 5))
N = T - (K*p) - 1 - (K*h) - (K - K*h +1)/2
df1 = h * K^2
df2 = (N * s) - ((K^2 * h) / 2) + 1
Qss = ((det(ΣUml) / det(ΣEml))^(1 / s) - 1) * (df2 / df1)
pvalss = ccdf(FDist(df1, df2), Qss)

# Test stat from Killian and Lütkepohl 2017, p.54
Q = T * (K - tr(ΣUml \ ΣEml))
df = h * K^2
pval = ccdf(Chisq(df), Q)
