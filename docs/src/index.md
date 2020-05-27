# VectorAutoregressions.jl

*VectorAutoregressions.jl* is a package for estimating Vector Autoregressions (VARs) using Julia.  The initial goal is to provide functionality comparable to the [`vars`](https://cran.r-project.org/package=vars) package in R.

Currently implemented:
- Calculate lag selection criteria.
- Estimate an unrestricted VAR using OLS.
- Check the stability of an estimated VAR.
