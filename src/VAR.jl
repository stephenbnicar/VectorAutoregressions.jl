"""
    VarEstimate

`struct` holding the output from a call to [`VAR`](@ref).
"""
struct VarEstimate
    "`DataFrame` or `TimeArray` of observations on endogenous variables"
    Y::Union{DataFrame,TimeArray}
    X::Union{DataFrame,TimeArray,Nothing}
    ynames::Array{String}
    xnames::Array{String}
    lags::Int
    constant::Bool
    trend::Bool
    obs::Int
    Z::Matrix
    B::Matrix
    seB::Matrix
    U::Matrix
    ΣU::Matrix
    Yhat::Matrix
end

"""
    VAR(data, lags; constant = true, trend = false) -> VarEstimate

Estimate an unrestricted vector autoregression (VAR) using OLS.

# Arguments
- `data` : `DataFrame` or `TimeArray` of observations on endogenous variables
- `lags::Int` : the number of lags
- `constant::Bool = true` : include an intercept term
- `trend::Bool = false` : include a linear trend
---
"""
function VAR(data::DataFrame, lags; constant::Bool = true, trend::Bool = false)
    if lags > size(data, 1)
        error("Number of lags is greater than number of observations.")
    end
    ynames = String.(names(data))
    datamat = Matrix(data)
    obs = size(datamat, 1) - lags
    Z, B, U, seB, ΣU, Yhat = varols(datamat, lags, constant, trend)
    VarEstimate(data, nothing, ynames, [""], lags, constant, trend,
        obs, Z, B, seB, U, ΣU, Yhat)
end

function VAR(data::TimeArray, lags; constant::Bool = true, trend::Bool = false)
    data = DataFrame(data)[:, 2:end]
    VAR(data, lags; constant = constant, trend = trend)
end

function varols(y, ylag, constant, trend)
    # Set up RHS Matrix
    Z = rhs_matrix(y, ylag, constant, trend)
    # Set up LHS Matrix
    Y = y[ylag+1:end, :]
    # Coefficient estimates
    B = (Z' * Z) \ (Z' * Y) # each column corresponds to an equation
    Yhat = Z * B
    U = Y - Yhat
    ΣU = (U' * U) / (size(Y, 1) - size(B, 1))
    ΣB = kron(ΣU, inv(Z' * Z))  # see Lutkepohl p.80
    seB = sqrt.(reshape(diag(ΣB), size(B, 1), size(B, 2)))
    return Z, B, U, seB, ΣU, Yhat
end

function show(io::IO, v::VarEstimate)
    ct = coeftable(v)
    K = length(ct)
    println(io, "VAR Estimation Results")
    println(io, "======================")
    print(io, "Endogenous variables: ")
    for k = 1:K-1
        print(io, "$(v.ynames[k]), ")
    end
    println(io, "$(v.ynames[K])")
    print(io, "Deterministic variables: ")
    v.constant && v.trend ? println(io, "constant, trend") :
    (v.constant ? println(io, "constant") : println(io))
    println(io, "Lags: $(v.lags)")
    println(io, "Sample size: $(v.obs)")
    for k = 1:K
        println(io)
        println(io, "Estimates for equation $(v.ynames[k]):")
        show(io, ct[k])
        println(io)
    end
end
