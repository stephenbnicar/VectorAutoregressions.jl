#==
This file uses the example introduced in Section 3.2.3 of Lutkepohl (2006)
Make sure to activate the VectorAutoregressions environment before running # hide
==#
using VectorAutoregressions
using CSV, DataFrames, DataFramesMeta
using Dates

exdir = dirname(dirname(pathof(VectorAutoregressions)))*"/examples";
datadf = CSV.read(exdir*"/lutkepohl_data.csv")

# Take the subset from 1960q1 to 1978q4, and use log-difference for each series
Y = @linq datadf |>
    where(:date .< Date(1979,3,1)) |>
    select(investment = diff(log.(:invest)), income = diff(log.(:income)),
        consumption = diff(log.(:cons)))

# ## Get Lag Selection Criteria
ls = lagselect(Y, 8)
# ## Estimate the VAR
v  = VAR(Y, ls.selection["AIC"])
# ## Check Stability of the VAR
sc = checkstable(v)
# ## Check for autocorrelation in the residuals
pt = portmanteau_test(v, 12)
