# Define structs for VAR models
mutable struct VarModel
    data
    lags::Int
    constant::Bool
    trend::Bool
end
