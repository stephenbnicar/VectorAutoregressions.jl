"""
A `struct` for estimating and storing a Vector Autoregression (VAR) estimated using ordinary least squares (OLS).

Constructors:

    VAR(Y, p; constant::Bool = true, trend::String = "none", exog = nothing, yvarnames = nothing)

Arguments:
---
* `Y` : Matrix or DataFrame of endogenous variables
* `p` : the number of lags
* `constant` : include a constant in the regression? (default is `true`)
* `trend` : "none", "linear", or "quadratic", as applicable (default is "none")
* `exog` : optional Vector, Matrix, or DataFrame of exogenous variables
* `yvarnames` : optional Vector of names for the endogenous variables

Fields (in addition to the arguments above):
---
* `Z` : Matrix of RHS variables used in the regression
* `B` : Matrix of coefficient estimates
* `A` : Array of coefficient estimates for lagged endogenous variables
* `U` : Matrix of residuals
* `SigmaU` : Residual variance/covariance matrix
* `seB` : Matrix of standard errors of coefficient estimates
---
"""
struct VAR <: AbstractVARModel
    Y::Matrix
    p::Int
    constant::Bool
    trend::String
    exog::Array{Float64}
    yvarnames::Vector{String}
    xvarnames::Vector{String}
    Z::Matrix
    B::Matrix
    A::Array{Float64,3}
    U::Matrix
    SigmaU::Matrix
    seB::Matrix

    function VAR(Y, p; constant::Bool = true, trend::String = "none", exog = nothing, yvarnames = nothing)
        if size(Y, 2) < 2
            error("A VAR requires at least two endogenous variables.")
        end
        if p > size(Y, 1)
            error("Number of lags is greater than number of observations.")
        end
        if trend ∉ ["none", "linear", "quadratic"]
            error("Invalid trend: must be 'none', 'linear', or 'quadratic'")
        end
        if isa(exog, Void)
            exog = Array{Float64}[]
        elseif size(exog, 1) !== size(Y, 1)
            error("Unequal number of observations for exogenous and endogenous variables.")
        end
        if isa(Y, DataFrame)
            yvarnames = string.(names(Y))
            Y = Matrix(Y)
        elseif isa(yvarnames, Void)
            yvarnames = [string("Y$i") for i in 1:size(Y,2)]
        end
        if isa(exog, DataFrame)
            xvarnames = string.(names(exog))
            exog = Matrix(exog)
        elseif size(exog, 1) == 0
            xvarnames = []
        else
            xvarnames = [string("X$i") for i in 1:size(exog,2)]
        end
        Z, B, A, U, SigmaU, seB = var_ols(Y, p, constant, trend, exog)
        new(Y, p, constant, trend, exog, yvarnames, xvarnames, Z, B, A, U, SigmaU, seB)
    end
end

"""
    var_ols(Y, p, constant, trend, exog)

Estimate a VAR using OLS. Used as part of the inner constructor for the [`VAR`](@ref) struct.
"""
function var_ols(Y, p, constant, trend, exog)
    nobs, K = size(Y)
    T = nobs - p
    Z = rhs_matrix(Y, p, constant, trend, exog)
    nparam = size(Z, 2)
    Y = Y[p+1:end, :]
    B = (Z'*Z)\(Z'*Y) # each column corresponds to an equation
    A = reshape(B[(nparam-K*p)+1:end, :]', K, K, p)
    U = Y - Z*B
    SigmaU = (U'*U)/(T - nparam)
    SigmaB = kron(SigmaU, inv(Z'*Z))  # following Lutkepohl(2005) p.80
    seB = sqrt.(reshape(diag(SigmaB), nparam, K))
    return Z, B, A, U, SigmaU, seB
end
