module VectorAutoregressions

using DataFrames, TimeSeries
using Distributions
using Dates, LinearAlgebra

import Base: show
import StatsBase: loglikelihood, aic, sample, coef, stderror, residuals

abstract type AbstractVarModel end

export  VarEstimate, VAR, IRF,
        LagSelectionCriteria,
        loglikelihood,
        aic, sic, hqc,
        coef, stderror, residuals,
        lagselect

include("matrix_utilities.jl")
include("VAR.jl")
include("IRF.jl")
include("diagnostics.jl")
include("lag_selection.jl")

end # module
