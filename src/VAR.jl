"""
    VAR(data, lags, constant::Bool = true, trend::Bool = false) -> VAR

Estimate an unrestricted vector autoregression (VAR) using OLS.

Arguments:
---
* `data` : `Matrix`, `DataFrame`, or `TimeArray` of observations on endogenous variables
* `lags` : the number of lags
* `constant` : boolean to indicate inclusion of intercept term (default is `true`)
* `trend` : boolean to indicate inclusion of a linear trend

Fields (in addition to the arguments above):
---
* `coef` : Matrix of coefficient estimates; each column corresponds to an equation
* `residuals` : Matrix of residuals; each column corresponds to an equation
* `se_coef` : Matrix of standard errors of coefficient estimates
* `vcov_residuals` : Variance-covariance Matrix{Float64} of residuals
---
"""
struct VAR <: AbstractVarModel
    data::Matrix{Float64}
    lags::Int
    constant::Bool
    trend::Bool
    coef::Matrix{Float64}
    se_coef::Matrix{Float64}
    residuals::Matrix{Float64}
    vcov_resid::Matrix{Float64}
end

function VAR(data::Matrix{Float64}, lags; constant::Bool = true, trend::Bool = false)
    if lags > size(data, 1)
        error("Number of lags is greater than number of observations.")
    end
    B, U, seB, Σᵤ = varols(data, lags, constant, trend)
    VAR(data, lags, constant, trend, B, seB, U, Σᵤ)
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
    # Set up RHS Matrix{Float64}
    Z = rhs_matrix(y, ylag, constant, trend)
    # Set up LHS Matrix{Float64}
    Y = y[ylag+1:end, :]
    # Coefficient estimates
    B = (Z' * Z) \ (Z' * Y) # each column corresponds to an equation
    U = Y - Z * B
    Σᵤ = (U' * U) / (size(Y, 1) - size(B, 1))
    ΣB = kron(Σᵤ, inv(Z' * Z))  # see Lutkepohl p.80
    seB = sqrt.(reshape(diag(ΣB), size(B, 1), size(B, 2)))
    return B, U, seB, Σᵤ
end
