#==

This file uses the example introduced in Section 3.2.3 of Lutkepohl (2006)

==#
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

ls = lagselect(Y, 8)
v  = VAR(Y, ls.selection[:AIC])
sc = checkstable(v)
