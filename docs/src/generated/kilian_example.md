```@meta
EditURL = "<unknown>/../examples/kilian_example.jl"
```

```@repl kilian_example
# Make sure to activate the VectorAutoregressions environment before running # hide
using CSV, DataFrames
using VectorAutoregressions

exdir = dirname(dirname(pathof(VectorAutoregressions)))*"/examples";
datadf = CSV.read(exdir*"/kilian_data.csv");
first(datadf, 6)
```

## Get Lag Selection Criteria

```@repl kilian_example
ls = lagselect(datadf, 8)
```

## Estimate the VAR

```@repl kilian_example
v  = VAR(datadf, ls.selection["HQC"])
```

## Check Stability of the VAR

```@repl kilian_example
vstable = checkstable(v)
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*
