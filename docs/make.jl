using Documenter, VectorAutoregressions

makedocs(
    modules=[VectorAutoregressions],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    sitename="VectorAutoregressions.jl",
    authors="Stephen Nicar",
)

deploydocs(
    repo="github.com/stephenbnicar/VectorAutoregressions.jl",
    # deploy_config = Documenter.GitHubActions()
)
