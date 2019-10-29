"""
    VAR(data, lags, constant::Bool = true, trend::Bool = false)

Estimate a vector autoregression (VAR) using OLS.

Arguments:
---
* `data` : Matrix or DataFrame of observations on endogenous variables
* `lags` : the number of lags
* `constant` : boolean to indicate inclusion of intercept term (default is `true`)
* `trend` : boolean to indicate inclusion of a linear trend

Fields (in addition to the arguments above):
---
* `coef` : Matrix of coefficient estimates; each column corresponds to an equation
* `residuals` : Matrix of residuals; each column corresponds to an equation
* `se_coef` : Matrix of standard errors of coefficient estimates
* `vcov_residuals` : Variance-covariance matrix of residuals
---
"""
struct VAR <: AbstractVarModel
    data::Matrix
    lags::Int
    constant::Bool
    trend::Bool
    coef::Matrix
    se_coef::Matrix
    residuals::Matrix
    vcov_resid::Matrix
end

function VAR(data::Matrix, lags; constant::Bool = true, trend::Bool = false)
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
    # Set up RHS Matrix
    Z = rhs_matrix(y, ylag, constant, trend)
    # Set up LHS Matrix
    Y = y[ylag+1:end, :]
    # Coefficient estimates
    B = (Z' * Z) \ (Z' * Y) # each column corresponds to an equation
    U = Y - Z * B
    Σᵤ = (U' * U) / (size(Y, 1) - size(B, 1))
    ΣB = kron(Σᵤ, inv(Z' * Z))  # see Lutkepohl p.80
    seB = sqrt.(reshape(diag(ΣB), size(B, 1), size(B, 2)))
    return B, U, seB, Σᵤ
end
