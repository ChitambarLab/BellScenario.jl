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
            "Types" => "BellScenario/types.md",
            "Scenarios" => "BellScenario/scenarios.md",
            "Adjacency Decomposition" => "BellScenario/adjacency_decomposition.md",
        ]
    ],
    repo="https://github.com/ChitambarLab/BellScenario.jl/blob/{commit}{path}#L{line}",
    sitename="BellScenario.jl",
    authors="Brian Doolittle",
    assets=String[],
)

deploydocs(;
    repo="github.com/ChitambarLab/BellScenario.jl",
)
