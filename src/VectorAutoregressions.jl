module VectorAutoregressions

using DataFrames, TimeSeries
using Distributions: ccdf, TDist
using Dates, LinearAlgebra
using StatsBase: CoefTable, autocov, autocor, crosscov, crosscor

import Base: show
import StatsBase: loglikelihood, aic, sample, coef, stderror, residuals,
    coeftable, fitted

export VarEstimate,
    VAR,
    VarLagSelection,
    lagselect,
    VarStabilityCheck,
    checkstable,
    portmanteau_test,
    loglikelihood,
    aic,
    sic,
    hqc,
    coef,
    stderror,
    residuals,
    coeftable,
    fitted

include("matrix_utilities.jl")
include("VAR.jl")
include("statsbase.jl")
include("diagnostics.jl")
include("lag_selection.jl")

end # module
