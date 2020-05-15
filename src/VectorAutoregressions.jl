module VectorAutoregressions

using DataFrames, TimeSeries
using Distributions: ccdf, TDist
using Dates, LinearAlgebra
using StatsBase: CoefTable, autocov, autocor, crosscov, crosscor

import Base: show
import StatsBase: loglikelihood, aic, sample, coef, stderror, residuals, coeftable

export VarEstimate,
    VAR,
    LagSelectionCriteria,
    StabilityCheck,
    loglikelihood,
    aic,
    sic,
    hqc,
    coef,
    stderror,
    residuals,
    coeftable,
    lagselect,
    checkstable

include("matrix_utilities.jl")
include("VAR.jl")
include("statsbase.jl")
include("diagnostics.jl")
include("lag_selection.jl")

end # module
