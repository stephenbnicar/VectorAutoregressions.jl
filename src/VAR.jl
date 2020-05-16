struct VarEstimate
    data::Union{DataFrame,TimeArray}
    varnames::Array{String}
    lags::Int
    constant::Bool
    trend::Bool
    obs::Int
    B::Matrix
    seB::Matrix
    U::Matrix
    ΣU::Matrix
end

"""
    VAR(data, lags, constant::Bool = true, trend::Bool = false) -> VarEstimate

Estimate an unrestricted vector autoregression (VAR) using OLS.

Arguments:
---
* `data` : `DataFrame` or `TimeArray` of observations on endogenous variables
* `lags` : the number of lags
* `constant` : `Bool` to indicate inclusion of intercept term
* `trend` : `Bool` to indicate inclusion of a linear trend
---
"""
function VAR(data::DataFrame, lags; constant::Bool = true, trend::Bool = false)
    if lags > size(data, 1)
        error("Number of lags is greater than number of observations.")
    end
    varnames = String.(names(data))
    datamat = Matrix(data)
    obs = size(datamat, 1) - lags
    B, U, seB, ΣU = varols(datamat, lags, constant, trend)
    VarEstimate(data, varnames, lags, constant, trend, obs, B, seB, U, ΣU)
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
    U = Y - Z * B
    ΣU = (U' * U) / (size(Y, 1) - size(B, 1))
    ΣB = kron(ΣU, inv(Z' * Z))  # see Lutkepohl p.80
    seB = sqrt.(reshape(diag(ΣB), size(B, 1), size(B, 2)))
    return B, U, seB, ΣU
end

function show(io::IO, v::VarEstimate)
    ct = coeftable(v)
    K = length(ct)
    println(io, "VAR Estimation Results")
    println(io, "======================")
    print(io, "Endogenous variables: ")
    for k = 1:K-1
        print(io, "$(v.varnames[k]), ")
    end
    println(io, "$(v.varnames[K])")
    print(io, "Deterministic variables: ")
    v.constant && v.trend ? println(io, "constant, trend") :
    (v.constant ? println(io, "constant") : println(io))
    println(io, "Lags: $(v.lags)")
    println(io, "Sample size: $(v.obs)")
    for k = 1:K
        println(io)
        println(io, "Estimates for equation $(v.varnames[k]):")
        show(io, ct[k])
        println(io)
    end
end


## New VarEstimate:
# Y original data matrix of endog variables
# X original data matrix of exog variables
# Z rhs matrix
# varnames y
# varnames x
# lags::Int
# constant::Bool
# trend::Bool
# obs::Int
# B::Matrix
# seB::Matrix
# U::Matrix
# ΣU::Matrix
# predicted values
# A matrix? (coeffs without const, trend)
