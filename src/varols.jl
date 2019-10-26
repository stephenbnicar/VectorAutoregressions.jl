# Inputs:
# Matrix of endogenous vars
# Constant
# Trend
# Number of Lags
# Exogenous vars
# Exogenous lags

function varols(y::Matrix, ylag::Int; constant::Bool=true, trend::Bool=false)

    # Set up RHS Matrix
    rhs = rhs_matrix(y, ylag, constant, trend)

    # Set up LHS Matrix
    lhs = y[ylag+1:end, :]

    # Coefficient estimates
    B = (rhs'*rhs)\(rhs'*lhs) # each column corresponds to an equation
    A = B[(constant+trend+1):end, :]
    ν = constant == true ? B[1:constant+trend, :] : []
    U = lhs - rhs*B
    Σᵤ = (U'*U)/(size(lhs, 1) - size(B, 1))
    ΣB = kron(Σᵤ, inv(rhs'*rhs))  # see Lutkepohl p.80
    seB = sqrt.(reshape(diag(ΣB), size(B, 1), size(B, 2)))
    return ν, A, U, Σᵤ, seB
end
