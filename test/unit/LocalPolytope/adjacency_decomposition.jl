using Test

@testset "./src/LocalPolytope/adjacency_decomposition.jl" begin

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

end
