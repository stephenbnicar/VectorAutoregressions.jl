# Define structs for VAR models

abstract type VarModelEstimators end



mutable struct VarModel{E}

    data
    lags::Int
    constant::Bool
    trend::Bool
end
