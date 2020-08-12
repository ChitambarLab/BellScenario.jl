using Test, LinearAlgebra

@testset "./src/LocalPolytope/vertices.jl" begin

using BellScenario

@testset "vertices()" begin
    @testset "vertices(PrepareAndMeasure::Scenario)" begin
        @testset "simple cases" begin
            PM = PrepareAndMeasure(2,2,2)
            gen_vertices = LocalPolytope.vertices(PM, rep="generalized")
            norm_vertices = LocalPolytope.vertices(PM, rep="normalized")

            @test gen_vertices == [[1,0,1,0],[0,1,0,1],[1,0,0,1],[0,1,1,0]]
            @test map(v -> convert(DeterministicStrategy, v, PM, rep="generalized"), gen_vertices) == [
                [1 1;0 0],[0 0;1 1],[1 0;0 1],[0 1;1 0]
            ]

            @test norm_vertices == [[1,1],[0,0],[1,0],[0,1]]
            @test map(v -> convert(DeterministicStrategy, v, PM, rep="normalized"), norm_vertices) == [
                [1 1;0 0],[0 0;1 1],[1 0;0 1],[0 1;1 0]
            ]

            @test_throws DomainError LocalPolytope.vertices(PM, rep="no-signaling")
        end

        @testset "vertices compatible with conversion to DeterministicStrategy" begin
            PM = PrepareAndMeasure(5,5,5)
            vertices = LocalPolytope.vertices(PM)

            @test map(
                v -> convert(DeterministicStrategy, v, PM, rep="normalized"),
            vertices) isa Vector{DeterministicStrategy}
        end

        @testset "rank_d_only = true" begin
            PM = PrepareAndMeasure(4,4,3)
            vertices = LocalPolytope.vertices(PM, rank_d_only=true)
            det_strategies = map(v -> convert(DeterministicStrategy, v, PM, rep="normalized"), vertices)

            @test all(s -> rank(s) == 3, det_strategies)
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

        PM = PrepareAndMeasure(4,4,3)
        vertices = LocalPolytope.vertices(PM, rank_d_only=true)
        @test length(vertices) == LocalPolytope.num_vertices(PM, rank_d_only=true)
    end
end

end
