module VectorAutoregressions

using DataFrames, TimeSeries
using Distributions
using Dates, LinearAlgebra

import Base: show
import StatsBase: loglikelihood, aic, sample

abstract type AbstractVarModel end

export  VAR, IRF,
        LagSelectionCriteria,
        loglikelihood,
        aic, sic, hqc,
        lagselect

include("matrix_utilities.jl")
include("VAR.jl")
include("IRF.jl")
include("diagnostics.jl")
include("lag_selection.jl")

end # module
