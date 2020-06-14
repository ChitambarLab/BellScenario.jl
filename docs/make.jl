using Documenter, BellScenario

makedocs(;
    modules=[BellScenario],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/bdoolittle/BellScenario.jl/blob/{commit}{path}#L{line}",
    sitename="BellScenario.jl",
    authors="Brian Doolittle",
    assets=String[],
)

deploydocs(;
    repo="github.com/bdoolittle/BellScenario.jl",
)
