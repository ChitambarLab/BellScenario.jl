using Documenter, BellScenario

makedocs(;
    modules=[BellScenario],
    format=Documenter.HTML(
        assets  = ["assets/custom.css"],
    ),
    pages=[
        "Home" => "index.md",
        "User Guide" => "user_guide.md",
        "BellScenario.jl" => [
            "Overview" => "BellScenario/overview.md",
            "Scenarios" => "BellScenario/scenarios.md",
            "Strategies" => "BellScenario/strategies.md",
            "Quantum Strategies" => "BellScenario/quantum_strategies.md",
            "Games" => "BellScenario/games.md",
            "File I/O" => "BellScenario/file_io.md",
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
            "Quantum Optimization" => "Nonlocality/quantum_opt.md",
        ],
        "Development Manual" => "development_manual.md",
    ],
    repo="https://github.com/ChitambarLab/BellScenario.jl/blob/{commit}{path}#L{line}",
    sitename="BellScenario.jl",
    authors="Brian Doolittle",
    assets=String[],
)

deploydocs(;
    repo="github.com/ChitambarLab/BellScenario.jl",
)
