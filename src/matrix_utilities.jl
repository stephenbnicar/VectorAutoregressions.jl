

function lag_matrix(x::T, lag::Int) where {T<:AbstractMatrix}
    nobs, nx = size(x)
    xlag = zeros(nobs - lag, nx * lag)
    for i = 1:lag
        xlag[1:(nobs-lag), (nx*(i-1)+1):(nx*i)] = x[(lag+1-i):(nobs-i), 1:nx]
    end
    return xlag
end

function lag_matrix(x::T, lag::Int) where {T<:AbstractVector}
    x2 = reshape(x, length(x), 1)
    return lag_matrix(x2, lag)
end

"""
    companion(A, p) -> Matrix

Construct the companion matrix for VAR coefficient matrix `A`, where `A` has
    dimensions (Kp × K). K is the number of endogenous variables, p is the
    number of lags.
"""
function companion(A, p)
    K = size(A, 2)
    subI = [Matrix(1.0I, K * (p - 1), K * (p - 1)) zeros(K * (p - 1), K)]
    Acomp = vcat(A', subI)
    return Acomp
end

function rhs_matrix(y, ylag, constant, trend, x)
    rhs = lag_matrix(y, ylag)
    if !isa(x, Nothing)
        rhs = hcat(x[ylag+1:end, :], rhs)
    end
    if trend
        ltrend = collect(1.0:size(rhs, 1))
        rhs = hcat(ltrend, rhs)
    end
    if constant
        rhs = hcat(ones(Float64, size(rhs, 1)), rhs)
    end
    return rhs
end
