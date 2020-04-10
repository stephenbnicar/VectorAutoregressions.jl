cd(joinpath(@__DIR__, "../"))
using Pkg
pkg"activate ."
using CSV, DataFrames, TimeSeries
using VectorAutoregressions
using Dates

cd(@__DIR__)
data_df = CSV.read("kilian_data.csv")

# Lag Selection Criteria
ls_df = lagselect(data_df, 8)

# Estimate VAR
var_df = VAR(data_df[!, 1:3], ls_df.selection[:AIC])
ct = coeftable(var_df)

# Construct TimeArray from Data
# dates = collect(Date(1954, 12, 1):Dates.Month(3):Date(2007, 12, 1))
# data_df.timestamp = dates
# data_ts = TimeArray(data_df, timestamp = :timestamp)
# ls_ts = lagselect(data_ts, 8)
# var_ts = VAR(data_ts, ls_ts.selection[:AIC])
