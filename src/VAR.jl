struct VarEstimate
    data::Matrix
    lags::Int
    constant::Bool
    trend::Bool
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
* `data` : `Matrix`, `DataFrame`, or `TimeArray` of observations on endogenous variables
* `lags` : the number of lags
* `constant` : boolean to indicate inclusion of intercept term (default is `true`)
* `trend` : boolean to indicate inclusion of a linear trend
---
"""
function VAR(data::Matrix, lags; constant::Bool = true, trend::Bool = false)
    if lags > size(data, 1)
        error("Number of lags is greater than number of observations.")
    end
    B, U, seB, ΣU = varols(data, lags, constant, trend)
    VarEstimate(data, lags, constant, trend, B, seB, U, ΣU)
end

function VAR(data::DataFrame, lags; constant::Bool = true, trend::Bool = false)
    data = Matrix(data)
    VAR(data, lags; constant = constant, trend = trend)
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
    println(io, "VAR Estimation Results")
    println(io, "====================================")
    println(io, "Fields: data, lags, constant, trend,")
    println(io, "        B, seB, U, ΣU")
end




# varresult
# list of ‘lm’ objects.
#
# datamat
# The data matrix of the endogenous and explanatory variables.
#
# y
# The data matrix of the endogenous variables
#
# type
# A character, specifying the deterministic regressors.
#
# p
# An integer specifying the lag order.
#
# K
# An integer specifying the dimension of the VAR.
#
# obs
# An integer specifying the number of used observations.
#
# totobs
# An integer specifying the total number of observations.
#
# restrictions
# Either NULL or a matrix object containing the zero restrictions of the VAR(p).
#
# call
# The call to VAR().
