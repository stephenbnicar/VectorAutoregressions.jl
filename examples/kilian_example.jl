cd(joinpath(@__DIR__, "../"))
using Pkg
pkg"activate ."
using CSV, DataFrames, TimeSeries
using VectorAutoregressions
using Dates

cd(@__DIR__)
data_df = CSV.read("kilian_data.csv")
dates = collect(Date(1954,12,1):Dates.Month(3):Date(2007,12,1))
data_df.timestamp = dates

# TimeArray
data_ts = TimeArray(data_df, timestamp = :timestamp)

var_df = VAR(data_df[!, 1:3], 4)

var_ts = VAR(data_ts, 4)
