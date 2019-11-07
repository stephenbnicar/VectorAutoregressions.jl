module VectorAutoregressions

using DataFrames, TimeSeries
using Distributions
using Dates, LinearAlgebra

import Base: show
import StatsBase: loglikelihood, aic

abstract type AbstractVarModel end

export  VAR,
        LagSelectionCriteria,
        loglikelihood,
        aic, sic, hqc,
        lagselect

include("matrix_utilities.jl")
include("VAR.jl")
include("diagnostics.jl")
include("lag_selection.jl")

end # module
