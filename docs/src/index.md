# VectorAutoregressions.jl

The following example is from Section 3.2.3 of Lutkepohl (2006).

```julia
using VectorAutoregressions
using CSV
using DataFrames
using DataFramesMeta
```
The full dataset consists of West German fixed investment, disposable income, and consumption expenditures, quarterly and seasonally adjusted, for 1960 - 1982. For this example, only observations through 1978 Q4 are used. The data are in natural logs and are first-differenced.

```julia
datadf = CSV.read("lutkepohl_data.csv")
Y = @linq datadf |>
    where(:date .< Date(1979,3,1)) |>
    select(investment = diff(log.(:invest)), income = diff(log.(:income)),
        consumption = diff(log.(:cons)))
```
