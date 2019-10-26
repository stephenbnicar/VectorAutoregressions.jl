module VectorAutoregressions

using  DataFrames
using  Distributions, StatsBase
using  LinearAlgebra

abstract type AbstractVarModel end

export VAR

include("matrix_utilities.jl")
include("VAR.jl")

end # module
