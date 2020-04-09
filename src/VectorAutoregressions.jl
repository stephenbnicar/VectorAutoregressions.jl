module VectorAutoregressions

using DataFrames, TimeSeries
using Distributions
using Dates, LinearAlgebra
using StatsBase: CoefTable

import Base: show
import StatsBase: loglikelihood, aic, sample, coef, stderror, residuals, coeftable

abstract type AbstractVarModel end

export  VarEstimate, VAR, IRF,
        LagSelectionCriteria,
        loglikelihood,
        aic, sic, hqc,
        coef, stderror, residuals, coeftable,
        lagselect

include("matrix_utilities.jl")
include("VAR.jl")
include("statsbase.jl")
include("IRF.jl")
include("diagnostics.jl")
include("lag_selection.jl")

end # module
