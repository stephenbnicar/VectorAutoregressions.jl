module VectorAutoregressions

using DataFrames, TimeSeries
using Distributions: ccdf, TDist, Chisq, FDist
using Dates, LinearAlgebra
using StatsBase: CoefTable, autocov, autocor, crosscov, crosscor

import Base: show
import StatsBase: loglikelihood, aic, sample, coef, stderror, residuals, coeftable, fitted

export VAREstimate,
    VAR,
    LagSelection,
    StabilityCheck,
    PortmanteauTest,
    LMCorrTest,
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
include("IRF.jl")


end # module
