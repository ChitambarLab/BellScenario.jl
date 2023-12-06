using Test, XPORTA

@testset "./src/LocalPolytope/facets.jl" begin

using BellScenario

@testset "facets()" begin
    @testset "3-simplex normalized" begin
        vertices = [[1,0],[0,1],[0,0]]
        facet_dict = LocalPolytope.facets(vertices)

        @test facet_dict["equalities"] == Vector{Vector{Int64}}[]
        @test facet_dict["facets"] == [[-1,0,0],[0,-1,0],[1,1,1]]
        @test length(keys(facet_dict)) == 2
    end

    @testset "3-simplex generalized" begin
        vertices = [[1,0,0],[0,1,0],[0,0,1]]
        facet_dict = LocalPolytope.facets(vertices)

        @test facet_dict["equalities"] == [[1,1,1,1]]
        @test facet_dict["facets"] == [[0,-1,0,0],[0,0,-1,0],[0,1,1,1]]
        @test length(keys(facet_dict)) == 2
    end

    @testset "CH inequalities" begin
        scenario = BipartiteNonSignaling(2,2,2,2)
        vertices = LocalPolytope.vertices(scenario)
        facet_dict = LocalPolytope.facets(vertices)

        @test facet_dict["equalities"] == Vector{Vector{Int64}}[]
        @test facet_dict["facets"] == [
            [0, 0, 0, 0, -1, 0, 0, 0, 0], [0, 0, 0, 0, 0, -1, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, -1, 0, 0], [0, 0, 0, 0, 0, 0, 0, -1, 0],
            [-1, 0, 0, 0, 0, 1, 0, 0, 0], [-1, 0, 0, 0, 1, 0, 0, 0, 0],
            [0, -1, 0, 0, 0, 0, 0, 1, 0], [0, -1, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, -1, 0, 0, 0, 1, 0, 0], [0, 0, -1, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, -1, 0, 0, 0, 1, 0], [0, 0, 0, -1, 0, 1, 0, 0, 0],
            [-1, 0, -1, 0, 1, 1, 1, -1, 0], [-1, 0, 0, -1, 1, 1, -1, 1, 0],
            [0, -1, -1, 0, 1, -1, 1, 1, 0], [0, -1, 0, -1, -1, 1, 1, 1, 0],
            [0, 1, 0, 1, 0, 0, 0, -1, 1], [0, 1, 1, 0, 0, 0, -1, 0, 1],
            [1, 0, 0, 1, 0, -1, 0, 0, 1], [1, 0, 1, 0, -1, 0, 0, 0, 1],
            [0, 1, 0, 1, 1, -1, -1, -1, 1], [0, 1, 1, 0, -1, 1, -1, -1, 1],
            [1, 0, 0, 1, -1, -1, 1, -1, 1], [1, 0, 1, 0, -1, -1, -1, 1, 1]
        ]
        @test length(keys(facet_dict)) == 2
    end
end

@testset "facets(::Polyhedron)" begin
    scenario = LocalSignaling(3,3,2)
    local_poly = LocalPolytope.vrep(scenario)

    @test !hrepiscomputed(local_poly)

    facets = LocalPolytope.facets(local_poly)
    @test length(facets) == 15

    @test hrepiscomputed(local_poly)

    match_vertices = LocalPolytope.vertices(scenario)
    match_facets = LocalPolytope.facets(match_vertices)["facets"]

    @test facets == match_facets
end

end
