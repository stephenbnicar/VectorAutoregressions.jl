cd(@__DIR__)

using CSV, DataFrames, TimeSeries, GLM
using VectorAutoregressions
using Dates

data_df = CSV.read("kilian_data.csv")
data_mat = Matrix(data_df)

endog_vars = string.(names(data_df))

K = size(data_mat, 2)
lags = 1
lagvars = Vector{String}()
for l = 1:lags
    for k = 1:K
        push!(lagvars, string(endog_vars[k], ".l$l"))
    end
end

Z = VectorAutoregressions.lag_matrix(data_mat, lags)
rhs = DataFrame(Z)
names!(rhs, Symbol.(lagvars))

rhs_terms = term.(names(rhs))
lhs_term = term(names(data_df)[1])
fterm = lhs_term ~ term(1) + foldl(+, rhs_terms)

data_all = [data_mat[lags+1:end, :] Z]
data_all_df = DataFrame(data_all)
names!(data_all_df, [names(data_df); names(rhs)])

var_glm = lm(fterm, data_all_df)

var = VAR(data_df, 1)
