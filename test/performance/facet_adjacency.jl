using BellScenario

function profile_test()
    vertices = map( v -> convert.(Int64, v), filter(v -> !in(v, LocalPolytope.vertices((5,1),(1,5))), LocalPolytope.vertices((5,1),(1,5), dits=2)))
    scenario = LocalSignaling(5,5,2)
    println(length(vertices))

    BG = BellGame([1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1], 2)

    # skip = [BellGame([1 0 0 0 0;1 0 0 0 0;1 0 0 0 0;1 0 0 0 0;0 0 0 0 0], 1)]
    @time facets = LocalPolytope.adjacency_decomposition(vertices, BG, scenario)#,  skip_games=skip)
    println(length(facets))
end

using ProfileView
profile_test()  # run once to trigger compilation (ignore this one)
@profview profile_test()
