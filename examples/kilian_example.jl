cd(joinpath(@__DIR__, "../"))
using Pkg
pkg"activate ."
using CSV, DataFrames
using VectorAutoregressions
using Dates

cd(@__DIR__)
if "kilian_data.csv" ∉ readdir()
    include("process_kilian_data.jl")
end
datadf = CSV.read("kilian_data.csv")

# Lag Selection Criteria
ls = lagselect(datadf, 8)

# Estimate VAR
v  = VAR(datadf, ls.selection["AIC"])
# v  = VAR(datadf, 4)

# Check Stability
vstable = checkstable(v)
