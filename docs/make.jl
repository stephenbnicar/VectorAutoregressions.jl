using Documenter, VectorAutoregressions

makedocs(;
    modules=[VectorAutoregressions],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/stephenbnicar/VectorAutoregressions.jl/blob/{commit}{path}#L{line}",
    sitename="VectorAutoregressions.jl",
    authors="Stephen Nicar",
    assets=String[],
)

deploydocs(;
    repo="github.com/stephenbnicar/VectorAutoregressions.jl",
)
