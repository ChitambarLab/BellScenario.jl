using Test, XPORTA

using Polyhedra: HalfSpace, vrep, convexhull, points

@testset "./src/ConvexPolytope/adjacency_decomposition.jl" begin

using BellScenario: adjacent_facets, rotate_facet, adjacency_decomposition, LocalPolytope, PrepareAndMeasure, BellGame

@testset "rotate_facet()" begin
    @testset "simplex rotation" begin
        F = HalfSpace([-1,0,0], 0)
        G1 = HalfSpace([0,0,-1], 0)
        G2 = HalfSpace([0,0,1], 1)

        xbar = [1,0,0]

        # neighboring facet 1
        @test rotate_facet(F,G1,xbar) == HalfSpace([0,0,-1], 0)
        @test rotate_facet(F,G2,xbar) == HalfSpace([1,0,1], 1)
    end
end

@testset "adjacent_facets" begin
    @testset "41-2-14 polytope" begin
        vertices = map( v -> convert.(Int64, v), LocalPolytope.vertices((4,1),(1,4), dits=2))
        PM = PrepareAndMeasure(4,4,2)

        F = HalfSpace([1,0,0,-1,0,1,0,-1,0,0,1,-1], 1)
        adj_facets = adjacent_facets(F, vertices, PM)

        println(adj_facets)
    end
end

@testset "adjacency_decomposition" begin
    @testset  "41-2-14 polytope" begin
        vertices = map( v -> convert.(Int64, v), LocalPolytope.vertices((4,1),(1,4), dits=2))
        PM = PrepareAndMeasure(4,4,2)

        BG = BellGame([1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1], 2)

        facets = adjacency_decomposition(BG, vertices, PM)

        println(facets)
    end

    # @testset "51-2-15 polytope" begin
    #     vertices = map( v -> convert.(Int64, v), LocalPolytope.vertices((5,1),(1,5), dits=2))
    #     PM = PrepareAndMeasure(5,5,2)
    #
    #     BG = BellGame([1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1], 2)
    #
    #     facets = adjacency_decomposition(BG, vertices, PM)
    #
    #     println(facets)
    # end

    # @testset "81-2-14 polytope" begin
    #     vertices = map( v -> convert.(Int64, v), LocalPolytope.vertices((8,1),(1,4), dits=2))
    #     println(length(vertices))
    #     PM = PrepareAndMeasure(8,4,2)
    #
    #     BG = BellGame([1 1 0 0 0 0 0 0;1 0 1 0 0 0 0 0;0 1 1 0 0 0 0 0;0 0 0 1 0 0 0 0], 3)
    #
    #     facets = adjacency_decomposition(BG, vertices, PM)
    #
    #     println(facets)
    # end

    @testset "51-3-15 polytope" begin
        vertices = map( v -> convert.(Int64, v), LocalPolytope.vertices((6,1),(1,6), dits=4))
        println(length(vertices))
        PM = PrepareAndMeasure(6,6,4)

        BG = BellGame([1 0 0 0 0 0;1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0], 4)

        facets = adjacency_decomposition(BG, vertices, PM)

        println(facets)
    end
end

end
