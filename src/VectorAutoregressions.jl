module VectorAutoregressions

using DataFrames, TimeSeries
using Distributions, StatsBase
using Dates, LinearAlgebra

abstract type AbstractVarModel end

export VAR

include("matrix_utilities.jl")
include("VAR.jl")

end # module
