"""
    VAREstimate

# Fields
- `Y::Union{DataFrame,TimeArray}`
- `X::Union{DataFrame,TimeArray,Nothing}`
- `ynames::Array{String}`
- `xnames::Array{String}`
- `lags::Int`
- `constant::Bool`
- `trend::Bool`
- `obs::Int`
- `Z::Matrix`
- `B::Matrix`
- `seB::Matrix`
- `U::Matrix`
- `ΣU::Matrix`
- `Yhat::Matrix`
"""
struct VAREstimate
    Y::Union{DataFrame,TimeArray}
    X::Union{DataFrame,TimeArray,Nothing}
    ynames::Array{String}
    xnames::Array{String}
    lags::Int
    constant::Bool
    trend::Bool
    obs::Int
    Z::Matrix
    B::Matrix
    seB::Matrix
    U::Matrix
    ΣU::Matrix
    Yhat::Matrix
end

"""
    VAR(endog, lags; constant = true, trend = false, exog = nothing) -> VAREstimate

Estimate an unrestricted vector autoregression (VAR) using OLS.

# Arguments
- `endog` : `DataFrame` or `TimeArray` of observations on endogenous variables
- `lags::Int` : the number of lags
- `constant::Bool = true` : include an intercept term
- `trend::Bool = false` : include a linear trend
- `exog` : `DataFrame` or `TimeArray` of observations on exogenous variables
"""
function VAR(
    endog::DataFrame,
    lags;
    constant::Bool = true,
    trend::Bool = false,
    exog::Union{DataFrame,Nothing} = nothing,
)
    if lags > size(endog, 1)
        error("Number of lags is greater than number of observations.")
    end
    if !isa(exog, Nothing) && size(exog, 1) !== size(endog, 1)
        error("Unequal number of observations for exogenous and endogenous variables.")
    end
    ynames = String.(names(endog))
    xnames = !isa(exog, Nothing) ? String.(names(exog)) : [""]
    datamat = Matrix(endog)
    xmat = !isa(exog, Nothing) ? Matrix(exog) : nothing
    obs = size(datamat, 1) - lags
    Z, B, U, seB, ΣU, Yhat = varols(datamat, lags, constant, trend, xmat)
    VAREstimate(
        endog,
        exog,
        ynames,
        xnames,
        lags,
        constant,
        trend,
        obs,
        Z,
        B,
        seB,
        U,
        ΣU,
        Yhat,
    )
end

function VAR(
    endog::TimeArray,
    lags;
    constant::Bool = true,
    trend::Bool = false,
    exog = nothing,
)
    endog = DataFrame(endog)[:, 2:end]
    if !isa(exog, Nothing)
        exog = DataFrame(exog)[:, 2:end]
    end
    VAR(endog, lags; constant = constant, trend = trend, exog = exog)
end

function varols(y, ylag, constant, trend, x)
    # Set up RHS Matrix
    Z = rhs_matrix(y, ylag, constant, trend, x)
    # Set up LHS Matrix
    Y = y[ylag+1:end, :]
    # Coefficient estimates: each column of B corresponds to an equation
    B = (Z' * Z) \ (Z' * Y)
    Yhat = Z * B
    U = Y - Yhat
    ΣU = (U' * U) / (size(Y, 1) - size(B, 1))
    ΣB = kron(ΣU, inv(Z' * Z))  # see Lutkepohl p.80
    seB = sqrt.(reshape(diag(ΣB), size(B, 1), size(B, 2)))
    return Z, B, U, seB, ΣU, Yhat
end

function show(io::IO, v::VAREstimate)
    ct = coeftable(v)
    K = length(ct)
    println(io, typeof(v))
    println(io)
    print(io, "Endogenous variables: ")
    for k = 1:K-1
        print(io, "$(v.ynames[k]), ")
    end
    println(io, "$(v.ynames[K])")
    print(io, "Deterministic variables: ")
    if !isa(v.X, Nothing)
        for xn = 1:length(v.xnames)
            print(io, "$(v.xnames[xn]), ")
        end
    end
    v.constant && v.trend ? println(io, "intercept, trend") :
    (v.constant ? println(io, "intercept") : println(io))
    println(io, "Lags: $(v.lags)")
    println(io, "Sample size: $(v.obs)")
    for k = 1:K
        println(io)
        println(io, "Estimates for equation $(v.ynames[k]):")
        show(io, ct[k])
        println(io)
    end
end
