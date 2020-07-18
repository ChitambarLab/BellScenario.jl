using BellScenario: adjacency_decomposition, LocalPolytope, PrepareAndMeasure, BellGame

function profile_test()
    vertices = map( v -> convert.(Int64, v), LocalPolytope.vertices((5,1),(1,5), dits=2))
    PM = PrepareAndMeasure(5,5,2)

    BG = BellGame([1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1], 2)

    @time adjacency_decomposition(BG, vertices, PM)

end

using ProfileView
@profview profile_test()  # run once to trigger compilation (ignore this one)
@profview profile_test()
