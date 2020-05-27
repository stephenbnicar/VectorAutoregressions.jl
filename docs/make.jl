using Documenter
using VectorAutoregressions
using Literate

EXAMPLE = dirname(dirname(pathof(VectorAutoregressions)))*"/examples/kilian_example.jl"
OUTPUT = joinpath(@__DIR__, "src/generated")
codefence = "```@repl kilian_example" => "```"
Literate.markdown(EXAMPLE, OUTPUT, codefence = codefence, documenter=true)

makedocs(
    modules=[VectorAutoregressions],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
        "Example" => "generated/kilian_example.md",
        "API"  => "api.md"
    ],
    sitename="VectorAutoregressions.jl",
    authors="Stephen Nicar",
)

deploydocs(
    repo="github.com/stephenbnicar/VectorAutoregressions.jl",
)
