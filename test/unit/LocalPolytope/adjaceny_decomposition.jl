using Test, XPORTA

@testset "./src/ConvexPolytope/adjacency_decomposition.jl" begin

using BellScenario

@testset "rotate_facet()" begin
    @testset "simplex rotation" begin
        F = [-1,0,0,0]
        G1 = [0,0,-1,0]
        G2 = [0,0,1,1]

        xbar = [1,0,0]

        # neighboring facets
        @test LocalPolytope.rotate_facet(F,G1,xbar) == [0,0,-1,0]
        @test LocalPolytope.rotate_facet(F,G2,xbar) == [1,0,1,1]
    end
end

@testset "adjacent_facets" begin
    @testset "41-2-14 polytope" begin
        PM = PrepareAndMeasure(4,4,2)
        vertices = LocalPolytope.vertices(PM)

        F = [1,0,0,0,1,0,0,0,1,-1,-1,-1,1]
        adj_facets = LocalPolytope.adjacent_facets(vertices, F)

        @test length(adj_facets) == 48

        adj_games = map(f -> convert(BellGame, f, PM, rep="normalized"), adj_facets)
        canonical_games = unique(map(g -> LocalPolytope.generator_facet(g, PM), adj_games))

        @test canonical_games == [
            [1 0 0 0;1 0 0 0;1 0 0 0;0 0 0 0],
            [1 0 0 0;1 0 0 0;0 1 0 0;0 0 1 0],
            [2 0 0 0;1 1 1 0;0 2 0 0;0 0 1 1],
            [1 1 0 0;1 0 1 0;0 1 1 0;0 0 0 1],
            [2 0 0 0;1 1 1 0;0 2 0 0;0 0 2 0]
        ]
    end
end

@testset "adjacency_decomposition" begin
    @testset  "41-2-14 polytope skip positivity" begin
        PM = PrepareAndMeasure(4,4,2)
        vertices = LocalPolytope.vertices(PM)

        BG = BellGame([1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1], 2)

        skip = [BellGame([1 0 0 0;1 0 0 0;1 0 0 0;0 0 0 0],1)]
        dict = LocalPolytope.adjacency_decomposition(vertices, BG, PM, skip_games=skip)

        # positvity included, but is skipped in computation
        @test collect(keys(dict)) == [
            [1 0 0 0;1 0 0 0;0 1 0 0;0 0 1 0],
            [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1],
            [1 0 0 0;1 0 0 0;1 0 0 0;0 0 0 0],
            [2 0 0 0;1 1 1 0;0 2 0 0;0 0 1 1],
            [2 0 0 0;1 1 1 0;0 2 0 0;0 0 2 0],
            [1 1 0 0;1 0 1 0;0 1 1 0;0 0 0 1],
        ]

        @test dict[[1 0 0 0;1 0 0 0;1 0 0 0;0 0 0 0]]["skipped"]
    end

    @testset  "41-2-14 polytope no skips" begin
        PM = PrepareAndMeasure(4,4,2)
        vertices = LocalPolytope.vertices(PM)

        BG = BellGame([1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1], 2)

        dict = LocalPolytope.adjacency_decomposition(vertices, BG, PM)

        @test collect(keys(dict)) == [
            [1 0 0 0;1 0 0 0;0 1 0 0;0 0 1 0],
            [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1],
            [1 0 0 0;1 0 0 0;1 0 0 0;0 0 0 0],
            [2 0 0 0;1 1 1 0;0 2 0 0;0 0 1 1],
            [2 0 0 0;1 1 1 0;0 2 0 0;0 0 2 0],
            [1 1 0 0;1 0 1 0;0 1 1 0;0 0 0 1]
        ]

        @test dict[[1 0 0 0;1 0 0 0;1 0 0 0;0 0 0 0]]["skipped"] == false
    end

    @testset  "41-2-14 polytope filter out positivity vertices" begin
        PM = PrepareAndMeasure(4,4,2)
        vertices = LocalPolytope.vertices(PM, rank_d_only = true)

        BG = BellGame([1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1], 2)

        dict = LocalPolytope.adjacency_decomposition(vertices, BG, PM)

        @test collect(keys(dict)) == [
            [1 0 0 0;1 0 0 0;0 1 0 0;0 0 1 0],
            [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1],
            [1 0 0 0;1 0 0 0;1 0 0 0;0 0 0 0],
            [2 0 0 0;1 1 1 0;0 2 0 0;0 0 1 1],
            [2 0 0 0;1 1 1 0;0 2 0 0;0 0 2 0],
            [1 1 0 0;1 0 1 0;0 1 1 0;0 0 0 1],
            [1 1 1 1;0 0 0 0;0 0 0 0;0 0 0 0],
        ]
    end
end

end
