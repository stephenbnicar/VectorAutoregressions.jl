## Make sure to activate the VectorAutoregressions environment before running # hide
using CSV, DataFrames
using VectorAutoregressions

exdir = pkgdir(VectorAutoregressions) * "/examples";
datadf = DataFrame!(CSV.File(exdir * "/kilian_data.csv"));
first(datadf, 6)

# ## Get Lag Selection Criteria
ls = LagSelection(datadf, 8)

# ## Estimate the VAR
v = VAR(datadf, ls.selection["AIC"])

# ## Check Stability of the VAR
vstable = StabilityCheck(v)
