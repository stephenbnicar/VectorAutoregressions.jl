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
end


function rhs_matrix(y, ylag, constant, trend)
    rhsy = lag_matrix(y, ylag)
    if trend
        ltrend = collect(1.0:size(Z, 1))
        rhsy = hcat(ltrend, rhsy)
    end
    if constant
        rhsy = hcat(ones(Float64, size(rhsy, 1)), rhsy)
    end
    return rhsy
end

function lag_matrix(x::T, lag::Int) where {T<:AbstractMatrix}
    nobs, nx = size(x)
    xlag = zeros(nobs - lag, nx*lag)
    for i = 1:lag
       xlag[1:(nobs - lag), (nx*(i - 1) + 1):(nx*i)] = x[(lag+1-i):(nobs-i), 1:nx]
    end
    return xlag
end

function lag_matrix(x::T, lag::Int) where {T<:AbstractVector}
    x2 = reshape(x, length(x), 1)
    return lag_matrix(x2, lag)
end
