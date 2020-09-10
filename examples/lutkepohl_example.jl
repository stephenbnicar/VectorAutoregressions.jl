# The following reproduces the example introduced in Section 3.2.3 of Lütkepohl (2006)

using VectorAutoregressions
using CSV, DataFrames, DataFramesMeta
using Dates

exdir = pkgdir(VectorAutoregressions) * "/examples";
datadf = DataFrame!(CSV.File(exdir * "/lutkepohl_data.csv"));
first(datadf, 6)

# Take the subset from 1960q1 to 1978q4, and use log-difference for each series
Y = @linq datadf |>
      where(:date .< Date(1979, 3, 1)) |>
      select(
          investment = diff(log.(:invest)),
          income = diff(log.(:income)),
          consumption = diff(log.(:cons)),
      );
first(Y, 6)

# ## Get lag selection criteria
ls = LagSelection(Y, 8)
# ## Estimate the VAR
v = VAR(Y, ls.selection["AIC"])
# ## Check stability of the VAR
sc = StabilityCheck(v)
# ## Check for autocorrelation in the residuals
pt = PortmanteauTest(v, 12)
