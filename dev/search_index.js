var documenterSearchIndex = {"docs":
[{"location":"#VectorAutoregressions.jl-1","page":"Home","title":"VectorAutoregressions.jl","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"VectorAutoregressions.jl is a package for estimating Vector Autoregressions (VARs) using Julia.  The initial goal is to provide functionality comparable to the vars package in R.","category":"page"},{"location":"#Estimating-a-VAR-1","page":"Home","title":"Estimating a VAR","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"VAR\nVarEstimate","category":"page"},{"location":"#VectorAutoregressions.VAR","page":"Home","title":"VectorAutoregressions.VAR","text":"VAR(data, lags; constant = true, trend = false) -> VarEstimate\n\nEstimate an unrestricted vector autoregression (VAR) using OLS.\n\nArguments\n\ndata : DataFrame or TimeArray of observations on endogenous variables\nlags::Int : the number of lags\nconstant::Bool = true : include an intercept term\ntrend::Bool = false : include a linear trend\n\n\n\n\n\n\n\n","category":"function"},{"location":"#VectorAutoregressions.VarEstimate","page":"Home","title":"VectorAutoregressions.VarEstimate","text":"VarEstimate\n\nstruct holding the output from a call to VAR.\n\n\n\n\n\n","category":"type"},{"location":"#Diagnostics-1","page":"Home","title":"Diagnostics","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"lagselect\nloglikelihood\naic\nsic\nhqc","category":"page"},{"location":"#VectorAutoregressions.lagselect","page":"Home","title":"VectorAutoregressions.lagselect","text":"lagselect(data, lags, constant::Bool = true,\n    trend::Bool = false) -> LagSelectionCriteria\n\nCalculate AIC, SIC, and HQC lag selection criteria for an unrestricted VAR.\n\nArguments:\n\ndata : DataFrame or TimeArray of observations on endogenous variables\nlags : maximum number of lags\nconstant : boolean to indicate inclusion of intercept term (default is true)\ntrend : boolean to indicate inclusion of a linear trend\n\n\n\n\n\n","category":"function"},{"location":"#StatsBase.loglikelihood","page":"Home","title":"StatsBase.loglikelihood","text":"loglikelihood(v::VarEstimate)\n\nReturn the log-likelihood for VAR model v.\n\n\n\n\n\n","category":"function"},{"location":"#StatsBase.aic","page":"Home","title":"StatsBase.aic","text":"aic(v::VarEstimate)\n\nReturn Akaike's Information Criterion for VAR model v.\n\n\n\n\n\n","category":"function"},{"location":"#VectorAutoregressions.sic","page":"Home","title":"VectorAutoregressions.sic","text":"sic(v::VarEstimate)\n\nReturn the Schwarz (Bayesian) Information Criterion for VAR model v.\n\n\n\n\n\n","category":"function"},{"location":"#VectorAutoregressions.hqc","page":"Home","title":"VectorAutoregressions.hqc","text":"hqc(v::VarEstimate)\n\nReturn the Hannan-Quinn Criterion for VAR model v.\n\n\n\n\n\n","category":"function"},{"location":"#Other-1","page":"Home","title":"Other","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"coef\nstderror\nresiduals\nfitted","category":"page"},{"location":"#StatsBase.coef","page":"Home","title":"StatsBase.coef","text":"coef(v::VarEstimate)\n\nReturn the matrix of coefficients for VAR model v.\n\n\n\n\n\n","category":"function"},{"location":"#StatsBase.stderror","page":"Home","title":"StatsBase.stderror","text":"stderror(v::VarEstimate)\n\nReturn the standard errors for the coefficients of VAR model v.\n\n\n\n\n\n","category":"function"},{"location":"#StatsBase.residuals","page":"Home","title":"StatsBase.residuals","text":"residuals(v::VarEstimate)\n\nReturn the matrix of residuals for VAR model v.\n\n\n\n\n\n","category":"function"},{"location":"#StatsBase.fitted","page":"Home","title":"StatsBase.fitted","text":"fitted(v::VarEstimate)\n\nReturn the fitted values for VAR model v.\n\n\n\n\n\n","category":"function"}]
}
