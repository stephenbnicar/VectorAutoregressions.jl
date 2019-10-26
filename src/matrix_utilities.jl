

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

"""
    comp_matrix(A, p)

Construct the companion matrix for `A`.
"""
function comp_matrix(A, p)
    K = size(A, 2)
    subI = [Matrix(1.0I, K*(p-1), K*(p-1)) zeros(K*(p-1), K)]
    Acomp = vcat(A', subI)
    return Acomp
end

function rhs_matrix(y, ylag, constant, trend)
    rhsy = lag_matrix(y, ylag)
    if trend
        ltrend = collect(1.0:size(rhsy, 1))
        rhsy = hcat(ltrend, rhsy)
    end
    if constant
        rhsy = hcat(ones(Float64, size(rhsy, 1)), rhsy)
    end
    return rhsy
end
