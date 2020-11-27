using Test, LinearAlgebra

@testset "./src/LocalPolytope/vertices.jl" begin

using BellScenario

@testset "vertices(scenario :: LocalSignaling)" begin
    @testset "vertices(LocalSignaling::Scenario)" begin
        @testset "simple cases" begin
            scenario = LocalSignaling(2,2,2)
            gen_vertices = LocalPolytope.vertices(scenario, rep="generalized")
            norm_vertices = LocalPolytope.vertices(scenario, rep="normalized")

            @test gen_vertices == [[1,0,1,0],[0,1,0,1],[1,0,0,1],[0,1,1,0]]
            @test map(v -> convert(DeterministicStrategy, v, scenario, rep="generalized"), gen_vertices) == [
                [1 1;0 0],[0 0;1 1],[1 0;0 1],[0 1;1 0]
            ]

            @test norm_vertices == [[1,1],[0,0],[1,0],[0,1]]
            @test map(v -> convert(DeterministicStrategy, v, scenario, rep="normalized"), norm_vertices) == [
                [1 1;0 0],[0 0;1 1],[1 0;0 1],[0 1;1 0]
            ]

            @test_throws DomainError LocalPolytope.vertices(scenario, rep="no-signaling")
        end

        @testset "vertices compatible with conversion to DeterministicStrategy" begin
            scenario = LocalSignaling(5,5,5)
            vertices = LocalPolytope.vertices(scenario)

            @test map(
                v -> convert(DeterministicStrategy, v, scenario, rep="normalized"),
            vertices) isa Vector{DeterministicStrategy}
        end

        @testset "rank_d_only = true" begin
            scenario = LocalSignaling(4,4,3)
            vertices = LocalPolytope.vertices(scenario, rank_d_only=true)
            det_strategies = map(v -> convert(DeterministicStrategy, v, scenario, rep="normalized"), vertices)

            @test all(s -> rank(s) == 3, det_strategies)
        end
    end
end

@testset "vertices(scenario :: BlackBox)" begin
    @testset "BlackBox(3,2)" begin
        gen_verts = LocalPolytope.vertices(BlackBox(3,2), rep="generalized")
        norm_verts = LocalPolytope.vertices(BlackBox(3,2))

        @test gen_verts == [
            [1,0,0,1,0,0],[0,1,0,1,0,0],[0,0,1,1,0,0],
            [1,0,0,0,1,0],[0,1,0,0,1,0],[0,0,1,0,1,0],
            [1,0,0,0,0,1],[0,1,0,0,0,1],[0,0,1,0,0,1]
        ]
        @test norm_verts == [
            [1,0,1,0],[0,1,1,0],[0,0,1,0],
            [1,0,0,1],[0,1,0,1],[0,0,0,1],
            [1,0,0,0],[0,1,0,0],[0,0,0,0]
        ]
    end

    @testset "DomainErrors" begin
        @test_throws DomainError LocalPolytope.vertices(BlackBox(2,2),rep="nonsignaling")
    end
end

@testset  "num_vertices()" begin
    @testset "num_vertices(LocalSignaling)" begin
        for X in 1:5
            for B in 1:5
                for d in 1:min(X,B)
                    scenario = LocalSignaling(X, B, d)
                    vertices = LocalPolytope.vertices(scenario)
                    @test length(vertices) == LocalPolytope.num_vertices(scenario)
                end
            end
        end

        scenario = LocalSignaling(4,4,3)
        vertices = LocalPolytope.vertices(scenario, rank_d_only=true)
        @test length(vertices) == LocalPolytope.num_vertices(scenario, rank_d_only=true)
    end
end

@testset "black_box_strategies()" begin
    @test LocalPolytope.black_box_strategies(BlackBox(1,5)) == [[1 1 1 1 1]]
    @test LocalPolytope.black_box_strategies(BlackBox(2,2)) == [
        [1 1;0 0],[0 1;1 0],[1 0;0 1],[0 0;1 1]
    ]
    @test LocalPolytope.black_box_strategies(BlackBox(2,3)) == [
        [1 1 1;0 0 0],[0 1 1;1 0 0],[1 0 1;0 1 0],[0 0 1;1 1 0],
        [1 1 0;0 0 1],[0 1 0;1 0 1],[1 0 0;0 1 1],[0 0 0;1 1 1]
    ]

    @testset "spot checks" begin
        strategies = LocalPolytope.black_box_strategies(BlackBox(5,7))

        @test length(strategies) == 5^7
        @test length(unique(strategies)) == 5^7
        @test all(s -> is_deterministic(s), strategies)
    end
end

end
