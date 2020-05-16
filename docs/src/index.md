# VectorAutoregressions.jl

*VectorAutoregressions.jl* is a package for estimating Vector Autoregressions (VARs) using Julia.  The initial goal is to provide functionality comparable to the [`vars`](https://cran.r-project.org/package=vars) package in R.


## Estimating a VAR

```@docs
VAR
VarEstimate
```

## Diagnostics

```@docs
lagselect
loglikelihood
aic
sic
hqc
```

## Other

```@docs
coef
stderror
residuals
fitted
```
