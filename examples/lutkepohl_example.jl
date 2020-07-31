# This file reproduces the example introduced in Section 3.2.3 of Lütkepohl (2006)

using VectorAutoregressions
using CSV, DataFrames, DataFramesMeta
using Dates

exdir = pkgdir(VectorAutoregressions) * "/examples";
datadf = DataFrame!(CSV.File(exdir * "/lutkepohl_data.csv"))

# Take the subset from 1960q1 to 1978q4, and use log-difference for each series
Y = @linq datadf |>
      where(:date .< Date(1979, 3, 1)) |>
      select(
          investment = diff(log.(:invest)),
          income = diff(log.(:income)),
          consumption = diff(log.(:cons)),
      )

# ## Get lag selection criteria
ls = lagselect(Y, 8)
# ## Estimate the VAR
v = VAR(Y, ls.selection["AIC"])
# ## Check stability of the VAR
sc = checkstable(v)
# ## Check for autocorrelation in the residuals
pt = portmanteau_test(v, 12)
bgt = bg_test(v, 12)
