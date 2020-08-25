# Activate VectorAutoregressions project first
cd(@__DIR__)
include("../examples/lutkepohl_example.jl")
using LinearAlgebra

v = VAR(Y,2)
Y = v.Y
p = v.lags
constant = v.constant
trend = v.trend
B = v.B
U = v.U
Σ = v.ΣU
K = size(B, 2)

nparam = size(B, 1)
A = B[(nparam-K*p)+1:end, :]
A = reshape(A', (K, K, p))
h = 3
#   phi(:,:,i) is all simple responses to all shocks at horizon i
phi = zeros(K, K, h + 1)
phi[:,:,1] = Matrix(1.0I, K, K)
phi[:,:,2] = phi[:,:,1] * A[:, :, 1]
for i = 3:h+1
    for j = 1:p
        phi[:, :, i] += phi[:, :, i-j] * A[:, :, j]
    end
end
