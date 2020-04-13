cd(joinpath(@__DIR__, "../"))
using Pkg
pkg"activate ."
using CSV, DataFrames, TimeSeries
using VectorAutoregressions
using Dates

cd(@__DIR__)
datadf = CSV.read("kilian_data.csv")

# Lag Selection Criteria
ls = lagselect(datadf, 8)

# Estimate VAR
v  = VAR(datadf, ls.selection[:AIC])

# Construct TimeArray from Data
# dates = collect(Date(1954, 12, 1):Dates.Month(3):Date(2007, 12, 1))
# datadf.timestamp = dates
# datats = TimeArray(datadf, timestamp = :timestamp)
# l2 = lagselect(datats, 8)
# v2 = VAR(datats, ls_ts.selection[:AIC])
