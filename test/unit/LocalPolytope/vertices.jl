using Test, LinearAlgebra

@testset "./src/LocalPolytope/vertices.jl" begin

using BellScenario

@testset "vertices()" begin
    @testset "vertices(PrepareAndMeasure)" begin
        @testset "simple cases" begin
            PM = PrepareAndMeasure(2,2,2)
            @test LocalPolytope.vertices(PM) == [[1 1;0 0],[0 0;1 1],[1 0;0 1],[0 1;1 0]]

            PM = PrepareAndMeasure(5,5,5)
            vertices = LocalPolytope.vertices(PM)
            @test all(is_deterministic.(vertices))
        end

        @testset "rank_d_only = true" begin
            PM = PrepareAndMeasure(4,4,3)
            vertices = LocalPolytope.vertices(PM, rank_d_only=true)
            @test all(v -> rank(v) == 3, vertices)
        end
    end
end

@testset  "num_vertices()" begin
    @testset "num_vertices(PrepareAndMeasure)" begin
        for X in 1:5
            for B in 1:5
                for d in 1:min(X,B)
                    PM = PrepareAndMeasure(X, B, d)
                    vertices = LocalPolytope.vertices(PM)
                    @test length(vertices) == LocalPolytope.num_vertices(PM)
                end
            end
        end
    end
end

end
