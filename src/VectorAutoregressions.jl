module VectorAutoregressions

using DataFrames, TimeSeries
using Distributions: ccdf, TDist, Chisq
using Dates, LinearAlgebra
using StatsBase: CoefTable, autocov, autocor, crosscov, crosscor

import Base: show
import StatsBase: loglikelihood, aic, sample, coef, stderror, residuals, coeftable, fitted

export VarEstimate,
    VAR,
    VarLagSelection,
    lagselect,
    VarStabilityCheck,
    checkstable,
    PortmanteauTest,
    portmanteau_test,
    BreuschGodfreyTest,
    bg_test,
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
