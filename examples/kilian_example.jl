cd(joinpath(@__DIR__, "../"))
using Pkg
pkg"activate ."
using VectorAutoregressions
using CSV
using DataFrames

cd(@__DIR__)
data = CSV.read("kilian_data.csv")

datamat = Matrix(data)

vout = varols(datamat, 4)
