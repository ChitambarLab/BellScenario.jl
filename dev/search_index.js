var documenterSearchIndex = {"docs":
[{"location":"BellScenario/overview/#Overview","page":"Overview","title":"Overview","text":"","category":"section"},{"location":"bell_comm_scenarios/#Bell-Communication-Scenarios","page":"Communication Scenarios","title":"Bell Communication Scenarios","text":"","category":"section"},{"location":"LocalPolytope/utils/","page":"Utilities","title":"Utilities","text":"CurrentModule = LocalPolytope","category":"page"},{"location":"LocalPolytope/utils/#Utilities-for-Facets-and-Vertices","page":"Utilities","title":"Utilities for Facets and Vertices","text":"","category":"section"},{"location":"LocalPolytope/utils/","page":"Utilities","title":"Utilities","text":"dimension","category":"page"},{"location":"LocalPolytope/utils/#BellScenario.LocalPolytope.dimension","page":"Utilities","title":"BellScenario.LocalPolytope.dimension","text":"dimension( vertices :: Vector{Vector{Int64}} ) :: Int64\n\nFor the provided vertices (points), finds the dimension of the affine space spanned by the vertices. This method computes the dimension of a polytope, facet, or collection of points.\n\nThis also accepts matrices and DeterministicStrategy types as arguments:\n\ndimension( vertices :: Vector{Matrix{Int64}} ) :: Int64\n\ndimension( vertices :: Vector{DeterministicStrategy} ) :: Int64\n\n\n\n\n\n","category":"function"},{"location":"BellScenario/scenarios/","page":"Scenarios","title":"Scenarios","text":"CurrentModule = BellScenario","category":"page"},{"location":"BellScenario/scenarios/#Black-box-Scenarios","page":"Scenarios","title":"Black-box Scenarios","text":"","category":"section"},{"location":"BellScenario/scenarios/","page":"Scenarios","title":"Scenarios","text":"Scenario\nBlackBox\nBipartite\nPrepareAndMeasure","category":"page"},{"location":"BellScenario/scenarios/#BellScenario.Scenario","page":"Scenarios","title":"BellScenario.Scenario","text":"An abstract type to represent general black-box scenarios.\n\n\n\n\n\n","category":"type"},{"location":"BellScenario/scenarios/#BellScenario.BlackBox","page":"Scenarios","title":"BellScenario.BlackBox","text":"BlackBox(num_in :: Int64, num_out :: Int64) <: Scenario\n\nA black-box scenario considering a single device. A black-box is an uncharacterized device with a finite number of classical inputs and outputs.\n\nA DomainError is throw if parameters num_in or num_out is less than 1.\n\n\n\n\n\n","category":"type"},{"location":"BellScenario/scenarios/#BellScenario.Bipartite","page":"Scenarios","title":"BellScenario.Bipartite","text":"Bipartite(\n    A :: Tuple{Int64, Int64},\n    B :: Tuple{Int64, Int64};\n    dits :: Int = 1,\n    bidirectional :: Bool = false\n) <: Scenario\n\nA black-box scenario with two devices and possible communication between the devices. The keyword parameter dits describes the number of dits used for communication and bidirectional describes whether communication is simultaneous between each party.\n\n\n\n\n\n","category":"type"},{"location":"BellScenario/scenarios/#BellScenario.PrepareAndMeasure","page":"Scenarios","title":"BellScenario.PrepareAndMeasure","text":"PrepareAndMeasure(\n    X :: Int,\n    B :: Int,\n    d :: Int,\n)\n\nA black-box scenario with X inputs, B outputs, and d-dits of communication from the input device to the output device.\n\n\n\n\n\n","category":"type"},{"location":"BellScenario/strategies/","page":"Strategies","title":"Strategies","text":"CurrentModule = BellScenario","category":"page"},{"location":"BellScenario/strategies/#Black-box-Strategies","page":"Strategies","title":"Black-box Strategies","text":"","category":"section"},{"location":"BellScenario/strategies/","page":"Strategies","title":"Strategies","text":"AbstractStrategy\nStrategy\nDeterministicStrategy\nis_deterministic\nstrategy_dims","category":"page"},{"location":"BellScenario/strategies/#BellScenario.AbstractStrategy","page":"Strategies","title":"BellScenario.AbstractStrategy","text":"A stochastic matrix which represents a map from input to output for a given Bell scenario.\n\nBase library extensions for  AbstractStrategy:\n\nTwo strategies can be multiplied together resulting in a new strategy, e.g. S1*S2 = S3.\n\n*(S1::AbstractStrategy, S2::AbstractStrategy) :: Strategy\n\nThe chained strategies can describe a new black-box scenario, therefore, scenario parameter can also be passed to the multiplication operator.\n\n*(S1::AbstractStrategy, S2::AbstractStrategy, scenario::Scenario) :: Strategy\n\n\n\n\n\n","category":"type"},{"location":"BellScenario/strategies/#BellScenario.Strategy","page":"Strategies","title":"BellScenario.Strategy","text":"A strategy matrix describing the statistical behavior of a black-box.\n\nStrategy(conditionals :: Matrix{<:Real}) <: AbstractMatrix{Float64}\n\nBy default, the constructor creates a strategy for a 'BlackBox' scenario. However, a Scenario can be passed to the Strategy constructor.\n\nStrategy(conditionals :: Matrix{<:Real}, scenario :: Scenario)\n\nA DomainError is thrown if the provided Scenario does not match the dimension of the conditionals matrix.\n\n\n\n\n\n","category":"type"},{"location":"BellScenario/strategies/#BellScenario.DeterministicStrategy","page":"Strategies","title":"BellScenario.DeterministicStrategy","text":"A strategy matrix describing the deterministic behavior of a black-box.\n\nDeterministicStrategy(conditionals :: Matrix{Int64}) <: AbstractMatrix{Int64}\n\nBy default, the constructor creates a strategy for a 'BlackBox' scenario. However, a Scenario can be passed to the DeterministicStrategy constructor.\n\nDeterministicStrategy(conditionals :: Matrix{Int}, scenario :: Scenario)\n\nA DomainError is thrown if the provided Scenario does not match the dimension of the conditionals matrix.\n\nA DomainError is thrown if the elements of conditionals are not 0 or 1.\n\nBase library Extensions for DeterministicStrategy:\n\nThe product of two deterministic strategies is a DeterministicStrategy.\n\n*(S1::DeterministicStrategy, S2::DeterministicStrategy) :: DeterministicStrategy\n\nDeterministic strategies correspond to vertices of the local polytope. A vertex is simply represented by a Vector{Int64} and the rep argument specifies whether the vertex is in the \"normalized\" or \"generalized\" representation.\n\nVertex (Vector{Int64}) -> DeterministicStrategy\n\nconvert(\n    ::Type{DeterministicStrategy},\n    vertex  :: Vector{Int64},\n    scenario :: Scenario;\n    rep = \"normalized\" :: String\n)\n\nDeterministicStrategy -> Vertex (Vector{Int64})\n\nconvert(\n    ::Type{Vector{Int64}},\n    strategy :: DeterministicStrategy;\n    rep = \"normalized\" :: String\n)\n\n\n\n\n\n","category":"type"},{"location":"BellScenario/strategies/#BellScenario.is_deterministic","page":"Strategies","title":"BellScenario.is_deterministic","text":"is_deterministic( strategy :: AbstractMatrix  ) :: Bool\n\nReturns true if all elements of strategy are either 0 or 1 and the matrix is a valid conditional probability distribution.\n\n\n\n\n\n","category":"function"},{"location":"BellScenario/strategies/#BellScenario.strategy_dims","page":"Strategies","title":"BellScenario.strategy_dims","text":"strategy_dims( scenario :: Scenario ) :: Tuple{Int, Int}\n\nReturns the dimensions of the Conditionals matrix describing a Strategy for the Scenario at hand. Each Scenario, can place unique constraints on the matrix dimensions, therfore, separate methods are called for each concrete Scenario.\n\n\n\n\n\n","category":"function"},{"location":"LocalPolytope/generators/","page":"Generators","title":"Generators","text":"CurrentModule = LocalPolytope","category":"page"},{"location":"LocalPolytope/generators/#Generating-Vertices-and-Facets","page":"Generators","title":"Generating Vertices and Facets","text":"","category":"section"},{"location":"LocalPolytope/generators/","page":"Generators","title":"Generators","text":"generator_facet","category":"page"},{"location":"LocalPolytope/generators/#BellScenario.LocalPolytope.generator_facet","page":"Generators","title":"BellScenario.LocalPolytope.generator_facet","text":"generator_facet( BG :: BellGame, PM :: PrepareAndMeasure ) :: BellGame\n\nFinds the generating facet for the provided BellGame. The generator is provided in lexicographic normal form. The generating facet is found recursively by an algorithm which sorts by lexicographic scores.\n\n\n\n\n\n","category":"function"},{"location":"quantum_nonlocality/#Quantum-Nonlocality","page":"Quantum Nonlocality","title":"Quantum Nonlocality","text":"","category":"section"},{"location":"LocalPolytope/vertices/","page":"Vertices","title":"Vertices","text":"CurrentModule = LocalPolytope","category":"page"},{"location":"LocalPolytope/vertices/#Vertex-Enumeration","page":"Vertices","title":"Vertex Enumeration","text":"","category":"section"},{"location":"LocalPolytope/vertices/","page":"Vertices","title":"Vertices","text":"vertices\nnum_vertices","category":"page"},{"location":"LocalPolytope/vertices/#BellScenario.LocalPolytope.vertices","page":"Vertices","title":"BellScenario.LocalPolytope.vertices","text":"vertices(\n    PM :: PrepareAndMeasure;\n    rep = \"normalized\" :: String\n    rank_d_only = false :: Bool\n) ::  Vector{Vector{Int64}}\n\nGenerates the deterministic strategies for the local polytope of PrepareAndMeasure scenarios.  The rank_d_only keyword arg specifies  whether to  exclude vertices which use fewer dits  of communication and thus have a  matrix rank less than d.\n\nwarning: Warning\nThe vertices computed in this method are vectorized directly from a strategy matrix by column-majorization. These vertices are distinct from those produced by older LocalPolytope.vertices() methods which are row-majorized. This functionality is a trial run of performance improvements to polytope computation.\n\n\n\n\n\n","category":"function"},{"location":"LocalPolytope/vertices/#BellScenario.LocalPolytope.num_vertices","page":"Vertices","title":"BellScenario.LocalPolytope.num_vertices","text":"num_vertices( PM :: PrepareAndMeasure; rank_d_only = false :: Bool ) :: Int64\n\nCounts the numbr of polytope vertices for the specified PrepareAndMeasure scenario. If rank_d_only = true, then only strategies using  d-dits are counted.\n\n\n\n\n\n","category":"function"},{"location":"development_manual/#Development-Manual","page":"Development Manual","title":"Development Manual","text":"","category":"section"},{"location":"development_manual/#Develop-BellComm.jl-Pkg","page":"Development Manual","title":"Develop BellComm.jl Pkg","text":"","category":"section"},{"location":"development_manual/","page":"Development Manual","title":"Development Manual","text":"This project is packaged using Pkg.jl. For code changes to be reflected in the packaged software, the package must be set to development mode. This is done by running the command:","category":"page"},{"location":"development_manual/","page":"Development Manual","title":"Development Manual","text":"julia -e 'using Pkg; Pkg.develop(\"BellComm\")'","category":"page"},{"location":"development_manual/#Run-Tests","page":"Development Manual","title":"Run Tests","text":"","category":"section"},{"location":"development_manual/","page":"Development Manual","title":"Development Manual","text":"Run the BellComm.jl package tests (continuous integration):","category":"page"},{"location":"development_manual/","page":"Development Manual","title":"Development Manual","text":"julia -e 'using Pkg;  Pkg.test(\"BellComm\")'","category":"page"},{"location":"development_manual/","page":"Development Manual","title":"Development Manual","text":"Download test dependencies into local environment:","category":"page"},{"location":"development_manual/","page":"Development Manual","title":"Development Manual","text":"julia --project=test/ -e 'using Pkg; Pkg.instantiate()'","category":"page"},{"location":"development_manual/","page":"Development Manual","title":"Development Manual","text":"Run tests from a dev environment:","category":"page"},{"location":"development_manual/","page":"Development Manual","title":"Development Manual","text":"julia test/path/to/test_file.jl","category":"page"},{"location":"development_manual/#Build-Docs","page":"Development Manual","title":"Build Docs","text":"","category":"section"},{"location":"development_manual/","page":"Development Manual","title":"Development Manual","text":"Deploy docs server on local machine:","category":"page"},{"location":"development_manual/","page":"Development Manual","title":"Development Manual","text":"cd docs/build/; python3 -m http.server --bind localhost","category":"page"},{"location":"development_manual/","page":"Development Manual","title":"Development Manual","text":"Build docs from BellComm.jl source (continuous integration):","category":"page"},{"location":"development_manual/","page":"Development Manual","title":"Development Manual","text":"Download docs dependencies into local environment:\n julia --project=docs/ -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'","category":"page"},{"location":"development_manual/","page":"Development Manual","title":"Development Manual","text":"Build docs:\n julia --project=docs/ docs/make.jl","category":"page"},{"location":"development_manual/","page":"Development Manual","title":"Development Manual","text":"Build docs from a dev environment:","category":"page"},{"location":"development_manual/","page":"Development Manual","title":"Development Manual","text":"julia -e 'using Pkg; Pkg.add(\"Documenter\")'","category":"page"},{"location":"development_manual/","page":"Development Manual","title":"Development Manual","text":"julia docs/make.jl","category":"page"},{"location":"BellScenario/file_io/","page":"File I/O","title":"File I/O","text":"CurrentModule = BellScenario","category":"page"},{"location":"BellScenario/file_io/#File-I/O","page":"File I/O","title":"File I/O","text":"","category":"section"},{"location":"BellScenario/file_io/","page":"File I/O","title":"File I/O","text":"pretty_print_txt","category":"page"},{"location":"BellScenario/file_io/#BellScenario.pretty_print_txt","page":"File I/O","title":"BellScenario.pretty_print_txt","text":"pretty_print_txt( bell_games :: Vector{BellGame}, filename :: String )\n\nPrints a set of BellGame's to filename.txt in a human-readable form.\n\n\n\n\n\n","category":"function"},{"location":"LocalPolytope/adjacency_decomposition/","page":"Adjacency Decomposition","title":"Adjacency Decomposition","text":"CurrentModule = LocalPolytope","category":"page"},{"location":"LocalPolytope/adjacency_decomposition/#Adjacency-Decomposition","page":"Adjacency Decomposition","title":"Adjacency Decomposition","text":"","category":"section"},{"location":"LocalPolytope/adjacency_decomposition/","page":"Adjacency Decomposition","title":"Adjacency Decomposition","text":"adjacency_decomposition\nadjacent_facets\nrotate_facet","category":"page"},{"location":"LocalPolytope/adjacency_decomposition/#BellScenario.LocalPolytope.adjacency_decomposition","page":"Adjacency Decomposition","title":"BellScenario.LocalPolytope.adjacency_decomposition","text":"adjacenecy_decomposition(\n    vertices :: Vector{Vector{Int64}},\n    BG_seed :: BellGame,\n    PM :: PrepareAndMeasure;\n    kwargs\n)\n\nGiven a polytpe represented by vertices, returns the complete set of canonical facets for prepare and measure scenario PM. The adjacencydecomposition algorithm requires a seeded vertex which is supplied with the `BGseed` argument. Facets are returned in the lexicographic normal form.\n\nReturns a dictionary where the keys are canonical BellGames and the value is a dictionary with keys\n\n\"considered\" => true, if the facet was considered in the algorithm.\n\"skipped\" => true, if the facet was skipped (not considered).\n\"num_vertices\" => number of vertices.\n\"norm_facet\" => facet vector representative (normalized rep) of canonical game\n\nKeyword  arguments kwargs\n\nskip_games = [] ::  Vector{BellGame} - List of games to skip.\nmax_vertices = 100 :: Int64 - The maximum number of vertices to allow in target facets.\ndir = \"./\" :: String- Directory in which to createporta_tmp/and.json` files.\nlog = false :: Bool - If true, the facet dictionary is  written to .json each iteration.\nlog_filename = \"adjacency_decomposition_now.json\" :: String\n\n\n\n\n\n","category":"function"},{"location":"LocalPolytope/adjacency_decomposition/#BellScenario.LocalPolytope.adjacent_facets","page":"Adjacency Decomposition","title":"BellScenario.LocalPolytope.adjacent_facets","text":"adjacent_facets(\n    vertices :: Vector{Vector{Int64}},\n    F :: Vector{int64};\n    dir = \"./\" :: String,\n    cleanup = true ::Bool\n) :: Vector{Vector{Int64}}\n\nFor the polytope represented by vertices, returns the set of facets adjacent to F.\n\nFacet vector F and the return facet vectors are assumed to ba in the normalized subspace.\n\nThe  dir argument specifies where to where to write files and directories from XPORTA.jl. If cleanup is true, then a  porta_tmp directory is created as a subdirectory of dir.\n\nIf cleanup is false, the created  porta_tmp directory is not removed.\n\n\n\n\n\n","category":"function"},{"location":"LocalPolytope/adjacency_decomposition/#BellScenario.LocalPolytope.rotate_facet","page":"Adjacency Decomposition","title":"BellScenario.LocalPolytope.rotate_facet","text":"rotate_facet(\n    F :: Vector{Int64},\n    G :: Vector{Int64},\n    xbar :: Vector{Int64}\n) :: Vector{Int64}\n\nPerforms a rotation of facet F relative to non-included vertex xbar and returns the rotated facet vector. F is a polytope facet, G is a subfacet of F and xbar is a vertex not contained by F. By construction, the rotated facet contains xbar.\n\n\n\n\n\n","category":"function"},{"location":"BellScenario/quantum_strategies/","page":"Quantum Strategies","title":"Quantum Strategies","text":"CurrentModule = BellScenario","category":"page"},{"location":"BellScenario/quantum_strategies/#Quantum-Strategies","page":"Quantum Strategies","title":"Quantum Strategies","text":"","category":"section"},{"location":"BellScenario/quantum_strategies/","page":"Quantum Strategies","title":"Quantum Strategies","text":"quantum_strategy","category":"page"},{"location":"BellScenario/quantum_strategies/#BellScenario.quantum_strategy","page":"Quantum Strategies","title":"BellScenario.quantum_strategy","text":"Constructs a strategy matrix given quantum states and measurements. The supported scenarios include:\n\nBlackBox scenarios\n\nquantum_strategy(\n    Π :: Observables.AbstractPOVM,\n    ρ_states :: Vector{<:States.AbstractDensityMatrix}\n) :: Strategy\n\nPrepareAndMeasure scenarios\n\nquantum_strategy(\n    Π :: Observables.AbstractPOVM,\n    ρ_states :: Vector{<:States.AbstractDensityMatrix},\n    PM :: PrepareAndMeasure\n) :: Strategy\n\nA DomainError is thrown if the provided states and measurements are not compatible with the specified scenario.\n\n\n\n\n\n","category":"function"},{"location":"user_guide/#User-Guide","page":"User Guide","title":"User Guide","text":"","category":"section"},{"location":"BellScenario/quantum_opt/","page":"Quantum Optimization","title":"Quantum Optimization","text":"CurrentModule = BellScenario","category":"page"},{"location":"BellScenario/quantum_opt/#Quantum-Optimization","page":"Quantum Optimization","title":"Quantum Optimization","text":"","category":"section"},{"location":"BellScenario/quantum_opt/","page":"Quantum Optimization","title":"Quantum Optimization","text":"optimize_measurement","category":"page"},{"location":"BellScenario/quantum_opt/#BellScenario.optimize_measurement","page":"Quantum Optimization","title":"BellScenario.optimize_measurement","text":"optimize_measurement(\n    game :: BellGame,\n    ρ_states :: Vector{<:States.AbstractDensityMatrix},\n    PM :: PrepareAndMeasure\n)\n\nPerform semi-definite programming to find the optimal quantum measurement which maximizes the score for the specified set of quantum  states ρ_states and BellGame.\n\n\n\n\n\n","category":"function"},{"location":"BellScenario/games/","page":"Games","title":"Games","text":"CurrentModule = BellScenario","category":"page"},{"location":"BellScenario/games/#Black-Box-Games","page":"Games","title":"Black-Box Games","text":"","category":"section"},{"location":"BellScenario/games/","page":"Games","title":"Games","text":"AbstractGame\nGame\nBellGame","category":"page"},{"location":"BellScenario/games/#BellScenario.AbstractGame","page":"Games","title":"BellScenario.AbstractGame","text":"Games represent a linear inequality.\n\n\n\n\n\n","category":"type"},{"location":"BellScenario/games/#BellScenario.Game","page":"Games","title":"BellScenario.Game","text":"Game(game::Matrix{T}, β::Real) <: AbstractMatrix{T}\n\nA Game is played by a Matrix and β is the bounding or winning score of the game. A Game is matrix representation of a linear inequality. Each element is a linear scale factor for an element of a strategy matrix.\n\nType parameter T is typically an Int64 or Float64.\n\n\n\n\n\n","category":"type"},{"location":"BellScenario/games/#BellScenario.BellGame","page":"Games","title":"BellScenario.BellGame","text":"BellGame(game::Matrix{Int64}, β::Int64, scenario::Scenario)\n\nA tight Bell inequality for the correlation polytope for scenario.\n\nBase library extensions for BellGame:\n\nConversions between BellGames and other polytope facet represenations. A facet vector is a column-major vectorization of a game matrix, with the bound placed in the last element of the vector.\n\nFacet (Vector{Int64}) -> BellGame\n\nconvert(\n    ::Type{BellGame},\n    facet::Vector{Int64},\n    scenario::Scenario;\n    rep = \"normalized\"::String\n)\n\nBellGame -> Facet (Vector{Int64})\n\nconvert(::Type{Vector{Int64}}, BG::BellGame; rep = \"normalized\" :: String)\n\nBellGame's to XPORTA.IEQ\n\nconvert(::Type{IEQ}, bell_games::Vector{BellGame}; rep = \"normalized\" :: String)\n\nXPORTA.IEQ to BellGame's\n\nconvert(\n    ::Type{Vector{BellGame}},\n    ieq::IEQ, scenario::Scenario;\n    rep = \"normalized\" :: String\n)\n\n\n\n\n\n","category":"type"},{"location":"local_bounds/#Local-Bounds","page":"Local Bounds","title":"Local Bounds","text":"","category":"section"},{"location":"#BellComm.jl","page":"Home","title":"BellComm.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Compute Bell inequalities and their quantum violations.","category":"page"},{"location":"","page":"Home","title":"Home","text":"BellComm.jl performs two tasks:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Computes the local polytope for classical bipartite Bell Scenarios enhanced with classical communication.\nOptimizes quantum states and measurements to find maximal quantum violations of the local bounds.","category":"page"}]
}