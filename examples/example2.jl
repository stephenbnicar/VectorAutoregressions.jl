cd(joinpath(@__DIR__, "../"))
using Pkg
pkg"activate ."

using VectorAutoregressions
using CSV, DataFrames, DataFramesMeta
using Dates

cd(@__DIR__)
datadf = CSV.read("lutkepohl_data.csv")

# Take the subset from 1960q1 to 1978q4, and use log-difference for each series
Y = @linq datadf |>
    where(:date .< Date(1979,3,1)) |>
    select(investment = diff(log.(:invest)), income = diff(log.(:income)),
        consumption = diff(log.(:cons)))

# Create quarterly dummies
q1 = repeat([1, 0, 0, 0], Int(size(datadf,1)/4))
q2 = repeat([0, 1, 0, 0], Int(size(datadf,1)/4))
q3 = repeat([0, 0, 1, 0], Int(size(datadf,1)/4))
q4 = repeat([0, 0, 0, 1], Int(size(datadf,1)/4))
qdums = hcat(q2, q3, q4) # Array
qdums_df = DataFrame(q2 = q2, q3 = q3, q4 = q4) # DataFrame
qdums_df = qdums_df[18:end, :]

lags = 2
constant = true
trend    = false

v = VAR(Y, lags, trend = trend, exog = qdums_df)
