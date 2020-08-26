using Test

@testset "./src/LocalPolytope/utils.jl" begin

using BellScenario

@testset "dimension()" begin

    @testset "vectors" begin
        vertices = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[1,1,0,0]]
        @test LocalPolytope.dimension(vertices) == 3
    end

    @testset "matrices" begin
        vertices = [[1 0;0 0],[0 0;0 1],[1 1;0 0],[1 1;0 1]]
        @test LocalPolytope.dimension(vertices) == 3
    end

    @testset "Deterministic Strategies" begin
        vertices = DeterministicStrategy.([[1 1;0 0],[0 0;1 1],[1 0;0 1]])
        @test LocalPolytope.dimension(vertices) == 2
    end
end

end
