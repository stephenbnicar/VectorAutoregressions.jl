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
* `constant` : boolean to indicate inclusion of intercept term (default is `true`)
* `trend` : boolean to indicate inclusion of a linear trend
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
    println(io, "VAR Estimation Results:")
    println(io, "=======================")
    print(io, "Endogenous variables: ")
    for k = 1:K
        print(io, "$(v.varnames[k]) ")
    end
    print(io, "\n")
    println(io, "Deterministic variables:")
    println(io, "Sample size: $(v.obs)")
    for k = 1:K
        println(io)
        println(io, "Estimates for equation $(v.varnames[k]):")
        show(io, ct[k])
        println(io)
    end
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

# VAR Estimation Results:
# =========================
# Endogenous variables: e, prod, rw, U
# Deterministic variables: const
# Sample size: 82
# Log Likelihood: -175.819
# Roots of the characteristic polynomial:
# 0.995 0.9081 0.9081 0.7381 0.7381 0.1856 0.1429 0.1429
# Call:
# VAR(y = Canada, p = 2, type = "const")
#
#
# Estimation results for equation e:
# ==================================
