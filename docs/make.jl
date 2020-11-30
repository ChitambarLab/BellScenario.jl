using Documenter, BellScenario

makedocs(;
    modules=[BellScenario],
    format=Documenter.HTML(
        assets  = ["assets/custom.css"],
    ),
    pages=[
        "Home" => "index.md",
        "Communication Scenarios" => "bell_comm_scenarios.md",
        "Development Manual" => "development_manual.md",
        "Local Bounds" => "local_bounds.md",
        "Quantum Nonlocality" => "quantum_nonlocality.md",
        "User Guide" => "user_guide.md",
        "BellScenario.jl" => [
            "Overview" => "BellScenario/overview.md",
            "Scenarios" => "BellScenario/scenarios.md",
            "Games" => "BellScenario/games.md",
            "Strategies" => "BellScenario/strategies.md",
            "Quantum Strategies" => "BellScenario/quantum_strategies.md",
            "File I/O" => "BellScenario/file_io.md",
        ],
        "LocalPolytope.jl" => [
            "Vertices"  => "LocalPolytope/vertices.md",
            "Facets" => "LocalPolytope/facets.md",
            "Generators" => "LocalPolytope/generators.md",
            "Adjacency Decomposition" => "LocalPolytope/adjacency_decomposition.md",
            "Utilities" => "LocalPolytope/utils.md",
        ],
    ],
    repo="https://github.com/ChitambarLab/BellScenario.jl/blob/{commit}{path}#L{line}",
    sitename="BellScenario.jl",
    authors="Brian Doolittle",
    assets=String[],
)

deploydocs(;
    repo="github.com/ChitambarLab/BellScenario.jl",
)
