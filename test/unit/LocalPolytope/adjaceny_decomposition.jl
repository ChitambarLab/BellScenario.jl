using Test, XPORTA

using Polyhedra: HalfSpace

@testset "./src/ConvexPolytope/adjacency_decomposition.jl" begin

using BellScenario: LocalPolytope, PrepareAndMeasure, BellGame

@testset "rotate_facet()" begin
    @testset "simplex rotation" begin
        F = HalfSpace([-1,0,0], 0)
        G1 = HalfSpace([0,0,-1], 0)
        G2 = HalfSpace([0,0,1], 1)

        xbar = [1,0,0]

        # neighboring facets
        @test LocalPolytope.rotate_facet(F,G1,xbar) == HalfSpace([0,0,-1], 0)
        @test LocalPolytope.rotate_facet(F,G2,xbar) == HalfSpace([1,0,1], 1)
    end
end

@testset "adjacent_facets" begin
    @testset "41-2-14 polytope" begin
        vertices = map( v -> convert.(Int64, v), LocalPolytope.vertices((4,1),(1,4), dits=2))
        PM = PrepareAndMeasure(4,4,2)

        F = HalfSpace([1,0,0,-1,0,1,0,-1,0,0,1,-1], 1)
        adj_facets = LocalPolytope.adjacent_facets(vertices, F, PM)

        @test unique(adj_facets) == [
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
        vertices = map( v -> convert.(Int64, v), LocalPolytope.vertices((4,1),(1,4), dits=2))
        PM = PrepareAndMeasure(4,4,2)

        BG = BellGame([1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1], 2)

        skip = [BellGame([1 0 0 0;1 0 0 0;1 0 0 0;0 0 0 0],1)]
        games = LocalPolytope.adjacency_decomposition(vertices, BG, PM, skip_games=skip)

        # positvity included, but is skipped in computation
        @test games == [
            [1 0 0 0;1 0 0 0;0 1 0 0;0 0 1 0],
            [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1],
            [1 0 0 0;1 0 0 0;1 0 0 0;0 0 0 0],
            [2 0 0 0;1 1 1 0;0 2 0 0;0 0 1 1],
            [2 0 0 0;1 1 1 0;0 2 0 0;0 0 2 0],
            [1 1 0 0;1 0 1 0;0 1 1 0;0 0 0 1],
        ]
    end

    @testset  "41-2-14 polytope no skips" begin
        vertices = map( v -> convert.(Int64, v), LocalPolytope.vertices((4,1),(1,4), dits=2))

        PM = PrepareAndMeasure(4,4,2)

        BG = BellGame([1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1], 2)

        games = LocalPolytope.adjacency_decomposition(vertices, BG, PM)

        @test games == [
            [1 0 0 0;1 0 0 0;0 1 0 0;0 0 1 0],
            [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1],
            [1 0 0 0;1 0 0 0;1 0 0 0;0 0 0 0],
            [2 0 0 0;1 1 1 0;0 2 0 0;0 0 1 1],
            [2 0 0 0;1 1 1 0;0 2 0 0;0 0 2 0],
            [1 1 0 0;1 0 1 0;0 1 1 0;0 0 0 1]
        ]
    end

    @testset  "41-2-14 polytope filter out positivity vertices" begin
        vertices = filter(v -> !in(v,[
                [1 1 1 1 0 0 0 0 0 0 0 0]',
                [0 0 0 0 1 1 1 1 0 0 0 0]',
                [0 0 0 0 0 0 0 0 1 1 1 1]',
                [0 0 0 0 0 0 0 0 0 0 0 0]'
            ])  , map( v -> convert.(Int64, v), LocalPolytope.vertices((4,1),(1,4), dits=2)))

        PM = PrepareAndMeasure(4,4,2)

        BG = BellGame([1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1], 2)

        games = LocalPolytope.adjacency_decomposition(vertices, BG, PM)

        @test games == [
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
