# VectorAutoregressions.jl

*Estimating Vector Autoregressions (VARs) using Julia.*

[![][docs-dev-img]][docs-dev-url]
[![][ci-img]][ci-url]

*VectorAutoregressions.jl* is an **in-development** package for estimating VARs using Julia.  The initial
goal is to provide functionality comparable to the [`vars`](https://cran.r-project.org/package=vars) package in R. Longer-term goals include adding Bayesian estimation.

Currently implemented:
- [x] Calculate lag selection criteria.
- [x] Estimate an unrestricted VAR using OLS.
- [x] Check the stability of an estimated VAR.
- [x] Portmanteau test for autocorellation in the residuals.
- [x] LM (Breusch-Godfrey) test for autocorellation in the residuals.

To do:
- [ ] Test for normally-distributed residuals.
- [ ] Calculate simple and orthogonalized impulse responses.
- [ ] Add support for linear restrictions.
- [ ] Add support for structural VARs.

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://stephenbnicar.github.io/VectorAutoregressions.jl/dev
[ci-img]: https://github.com/stephenbnicar/VectorAutoregressions.jl/workflows/CI/badge.svg
[ci-url]: https://github.com/stephenbnicar/VectorAutoregressions.jl/actions?workflow=CI
