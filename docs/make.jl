push!(LOAD_PATH,"../src/")

using Documenter, BellScenario

DocMeta.setdocmeta!(BellScenario, :DocTestSetup, :(using BellScenario); recursive=true)

makedocs(;
    modules=[BellScenario],
    format=Documenter.HTML(
        assets = ["assets/custom.css"],
    ),
    pages=[
        "Home" => "index.md",
        "User Guide" => "user_guide.md",
        "BellScenario.jl" => [
            "Overview" => "BellScenario/overview.md",
            "Scenarios" => "BellScenario/scenarios.md",
            "Strategies" => "BellScenario/strategies.md",
            "Games" => "BellScenario/games.md",
            "Combinatorics" => "BellScenario/combinatorics.md",
        ],
        "LocalPolytope.jl" => [
            "Overview" => "LocalPolytope/overview.md",
            "Vertices" => "LocalPolytope/vertices.md",
            "Facets" => "LocalPolytope/facets.md",
            "Generators" => "LocalPolytope/generators.md",
            "Adjacency Decomposition" => "LocalPolytope/adjacency_decomposition.md",
            "Utilities" => "LocalPolytope/utils.md",
        ],
        "Nonlocality.jl" => [
            "Overview" => "Nonlocality/overview.md",
            "Optimize Measurements" => "Nonlocality/optimize_measurements.md",
        ],
        "Development Manual" => "development_manual.md",
    ],
    repo="https://github.com/ChitambarLab/BellScenario.jl/blob/{commit}{path}#L{line}",
    sitename="BellScenario.jl",
    authors="Brian Doolittle",
    warnonly = true,
)

deploydocs(;
    repo="github.com/ChitambarLab/BellScenario.jl",
    devbranch = "main",
)
