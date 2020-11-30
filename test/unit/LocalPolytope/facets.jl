using Test

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
        scenario = BipartiteNoSignaling(2,2,2,2)
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

end
