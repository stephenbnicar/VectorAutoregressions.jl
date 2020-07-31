var documenterSearchIndex = {"docs":
[{"location":"generated/kilian_example/#","page":"-","title":"-","text":"EditURL = \"<unknown>/../examples/kilian_example.jl\"","category":"page"},{"location":"generated/kilian_example/#","page":"-","title":"-","text":"# Make sure to activate the VectorAutoregressions environment before running # hide\nusing CSV, DataFrames\nusing VectorAutoregressions\n\nexdir = dirname(dirname(pathof(VectorAutoregressions)))*\"/examples\";\ndatadf = CSV.read(exdir*\"/kilian_data.csv\");\nfirst(datadf, 6)","category":"page"},{"location":"generated/kilian_example/#Get-Lag-Selection-Criteria-1","page":"-","title":"Get Lag Selection Criteria","text":"","category":"section"},{"location":"generated/kilian_example/#","page":"-","title":"-","text":"ls = lagselect(datadf, 8)","category":"page"},{"location":"generated/kilian_example/#Estimate-the-VAR-1","page":"-","title":"Estimate the VAR","text":"","category":"section"},{"location":"generated/kilian_example/#","page":"-","title":"-","text":"v  = VAR(datadf, ls.selection[\"HQC\"])","category":"page"},{"location":"generated/kilian_example/#Check-Stability-of-the-VAR-1","page":"-","title":"Check Stability of the VAR","text":"","category":"section"},{"location":"generated/kilian_example/#","page":"-","title":"-","text":"vstable = checkstable(v)","category":"page"},{"location":"generated/kilian_example/#","page":"-","title":"-","text":"","category":"page"},{"location":"generated/kilian_example/#","page":"-","title":"-","text":"This page was generated using Literate.jl.","category":"page"},{"location":"api/#Estimating-a-VAR-1","page":"API","title":"Estimating a VAR","text":"","category":"section"},{"location":"api/#","page":"API","title":"API","text":"VAR\nVarEstimate","category":"page"},{"location":"api/#VectorAutoregressions.VAR","page":"API","title":"VectorAutoregressions.VAR","text":"VAR(endog, lags; constant = true, trend = false, exog = nothing) -> VarEstimate\n\nEstimate an unrestricted vector autoregression (VAR) using OLS.\n\nArguments\n\nendog : DataFrame or TimeArray of observations on endogenous variables\nlags::Int : the number of lags\nconstant::Bool = true : include an intercept term\ntrend::Bool = false : include a linear trend\nexog : DataFrame or TimeArray of observations on exogenous variables\n\n\n\n\n\n","category":"function"},{"location":"api/#VectorAutoregressions.VarEstimate","page":"API","title":"VectorAutoregressions.VarEstimate","text":"VarEstimate\n\nFields\n\nY::Union{DataFrame,TimeArray}\nX::Union{DataFrame,TimeArray,Nothing}\nynames::Array{String}\nxnames::Array{String}\nlags::Int\nconstant::Bool\ntrend::Bool\nobs::Int\nZ::Matrix\nB::Matrix\nseB::Matrix\nU::Matrix\nΣU::Matrix\nYhat::Matrix\n\n\n\n\n\n","category":"type"},{"location":"api/#Diagnostics-1","page":"API","title":"Diagnostics","text":"","category":"section"},{"location":"api/#","page":"API","title":"API","text":"lagselect\nVarLagSelection\ncheckstable\nVarStabilityCheck\nportmanteau_test\nPortmanteauTest\nloglikelihood\naic\nsic\nhqc","category":"page"},{"location":"api/#VectorAutoregressions.lagselect","page":"API","title":"VectorAutoregressions.lagselect","text":"lagselect(data, maxlag; constant = true, trend = false) -> VarLagSelection\n\nCalculate AIC, SIC, and HQC lag selection criteria for an unrestricted VAR.\n\nArguments\n\ndata : DataFrame or TimeArray of observations on endogenous variables\nmaxlag::Int : the maximum number of lags to estimate\nconstant::Bool = true : include an intercept term\ntrend::Bool = false : include a linear trend\n\n\n\n\n\n","category":"function"},{"location":"api/#VectorAutoregressions.VarLagSelection","page":"API","title":"VectorAutoregressions.VarLagSelection","text":"VarLagSelection\n\nFields\n\nmaxlag::Int\ntable::DataFrame\nselection::Dict\n\n\n\n\n\n","category":"type"},{"location":"api/#VectorAutoregressions.checkstable","page":"API","title":"VectorAutoregressions.checkstable","text":"checkstable(v::VarEstimate) -> VarStabilityCheck\n\nCheck the stability of the estimated VAR model v.\n\n\n\n\n\n","category":"function"},{"location":"api/#VectorAutoregressions.VarStabilityCheck","page":"API","title":"VectorAutoregressions.VarStabilityCheck","text":"VarStabilityCheck\n\nstruct to hold the results of checkstable.\n\nFields\n\nisstable::Bool\neigenvals::Array{Number}\neigenmod::Array{Float64}\n\n\n\n\n\n","category":"type"},{"location":"api/#VectorAutoregressions.portmanteau_test","page":"API","title":"VectorAutoregressions.portmanteau_test","text":"portmanteau_test(v::VarEstimate, h::Int) -> PortmanteauTest\n\nConduct a multivariate portmanteau test for residual autocorrelation up to lag     h for VAR model v. Implements the adjusted test described on p.171 of Lutkepohl (2006).\n\n\n\n\n\n","category":"function"},{"location":"api/#VectorAutoregressions.PortmanteauTest","page":"API","title":"VectorAutoregressions.PortmanteauTest","text":"PortmanteauTest\n\nstruct to hold the results of portmanteau_test.\n\nFields\n\nU::Matrix\nh::Int\nQ::Float64\ndf::Int\npval::Float64\n\n\n\n\n\n","category":"type"},{"location":"api/#StatsBase.loglikelihood","page":"API","title":"StatsBase.loglikelihood","text":"loglikelihood(v::VarEstimate)\n\nReturn the log-likelihood for VAR model v.\n\n\n\n\n\n","category":"function"},{"location":"api/#StatsBase.aic","page":"API","title":"StatsBase.aic","text":"aic(v::VarEstimate)\n\nReturn Akaike's Information Criterion for VAR model v.\n\n\n\n\n\n","category":"function"},{"location":"api/#VectorAutoregressions.sic","page":"API","title":"VectorAutoregressions.sic","text":"sic(v::VarEstimate)\n\nReturn the Schwarz (Bayesian) Information Criterion for VAR model v.\n\n\n\n\n\n","category":"function"},{"location":"api/#VectorAutoregressions.hqc","page":"API","title":"VectorAutoregressions.hqc","text":"hqc(v::VarEstimate)\n\nReturn the Hannan-Quinn Criterion for VAR model v.\n\n\n\n\n\n","category":"function"},{"location":"api/#Other-1","page":"API","title":"Other","text":"","category":"section"},{"location":"api/#","page":"API","title":"API","text":"coef\nstderror\nresiduals\nfitted","category":"page"},{"location":"api/#StatsBase.coef","page":"API","title":"StatsBase.coef","text":"coef(v::VarEstimate)\n\nReturn the matrix of coefficients for VAR model v.\n\n\n\n\n\n","category":"function"},{"location":"api/#StatsBase.stderror","page":"API","title":"StatsBase.stderror","text":"stderror(v::VarEstimate)\n\nReturn the standard errors for the coefficients of VAR model v.\n\n\n\n\n\n","category":"function"},{"location":"api/#StatsBase.residuals","page":"API","title":"StatsBase.residuals","text":"residuals(v::VarEstimate)\n\nReturn the matrix of residuals for VAR model v.\n\n\n\n\n\n","category":"function"},{"location":"api/#StatsBase.fitted","page":"API","title":"StatsBase.fitted","text":"fitted(v::VarEstimate)\n\nReturn the fitted values for VAR model v.\n\n\n\n\n\n","category":"function"},{"location":"generated/lutkepohl_example/#","page":"Example","title":"Example","text":"EditURL = \"https://github.com/stephenbnicar/VectorAutoregressions.jl/blob/master/examples/lutkepohl_example.jl\"","category":"page"},{"location":"generated/lutkepohl_example/#","page":"Example","title":"Example","text":"This file reproduces the example introduced in Section 3.2.3 of Lütkepohl (2006)","category":"page"},{"location":"generated/lutkepohl_example/#","page":"Example","title":"Example","text":"using VectorAutoregressions\nusing CSV, DataFrames, DataFramesMeta\nusing Dates\n\nexdir = pkgdir(VectorAutoregressions) * \"/examples\";\ndatadf = DataFrame!(CSV.File(exdir * \"/lutkepohl_data.csv\"))","category":"page"},{"location":"generated/lutkepohl_example/#","page":"Example","title":"Example","text":"Take the subset from 1960q1 to 1978q4, and use log-difference for each series","category":"page"},{"location":"generated/lutkepohl_example/#","page":"Example","title":"Example","text":"Y = @linq datadf |>\n      where(:date .< Date(1979, 3, 1)) |>\n      select(\n          investment = diff(log.(:invest)),\n          income = diff(log.(:income)),\n          consumption = diff(log.(:cons)),\n      )","category":"page"},{"location":"generated/lutkepohl_example/#Get-lag-selection-criteria-1","page":"Example","title":"Get lag selection criteria","text":"","category":"section"},{"location":"generated/lutkepohl_example/#","page":"Example","title":"Example","text":"ls = lagselect(Y, 8)","category":"page"},{"location":"generated/lutkepohl_example/#Estimate-the-VAR-1","page":"Example","title":"Estimate the VAR","text":"","category":"section"},{"location":"generated/lutkepohl_example/#","page":"Example","title":"Example","text":"v = VAR(Y, ls.selection[\"AIC\"])","category":"page"},{"location":"generated/lutkepohl_example/#Check-stability-of-the-VAR-1","page":"Example","title":"Check stability of the VAR","text":"","category":"section"},{"location":"generated/lutkepohl_example/#","page":"Example","title":"Example","text":"sc = checkstable(v)","category":"page"},{"location":"generated/lutkepohl_example/#Check-for-autocorrelation-in-the-residuals-1","page":"Example","title":"Check for autocorrelation in the residuals","text":"","category":"section"},{"location":"generated/lutkepohl_example/#","page":"Example","title":"Example","text":"pt = portmanteau_test(v, 12)\nbgt = bg_test(v, 12)","category":"page"},{"location":"generated/lutkepohl_example/#","page":"Example","title":"Example","text":"","category":"page"},{"location":"generated/lutkepohl_example/#","page":"Example","title":"Example","text":"This page was generated using Literate.jl.","category":"page"},{"location":"#VectorAutoregressions.jl-1","page":"Home","title":"VectorAutoregressions.jl","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"VectorAutoregressions.jl is a package for estimating Vector Autoregressions (VARs) using Julia.  The initial goal is to provide functionality comparable to the vars package in R.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Currently implemented:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Calculate lag selection criteria.\nEstimate an unrestricted VAR using OLS.\nCheck the stability of an estimated VAR.\nPortmanteau test for autocorellation in the residuals.","category":"page"}]
}
