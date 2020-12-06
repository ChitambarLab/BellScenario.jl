using Test

@testset "./src/deterministic_strategies.jl" begin

using BellScenario

@testset "DeterministicStrategy()" begin
    D = DeterministicStrategy([1.0 0.0;0.0 1.0;0 0])

    @test D isa DeterministicStrategy
    @test D.conditionals isa Matrix{Int64}
    @test D == [1 0;0 1;0 0]
    @test D.scenario == BlackBox(3,2)

    @test_throws DomainError DeterministicStrategy([0.5 0.5;0.5 0.5])
    @test_throws DomainError DeterministicStrategy([1 0;0 1], BlackBox(3,2))
end

@testset "is_deterministic()" begin
    @test is_deterministic([1 0;0 1])
    @test !is_deterministic([0.5 0.5;0.5 0.5])
    @test is_deterministic([1.0 1.0;0.0 0.0])
end

@testset "strategy conversions" begin
    @testset "DeterministicStrategy -> Vector{Int64}" begin
        s1 = DeterministicStrategy([1 1 1;0 0 0;0 0 0])
        s2 = DeterministicStrategy([1 0 0;0 1 0;0 0 1])

        norm_v1 = convert(Vector{Int64}, s1, rep="normalized")
        norm_v2 = convert(Vector{Int64}, s2, rep="normalized")
        gen_v = convert(Vector{Int64}, s1, rep="generalized")

        @test norm_v1 isa Vector{Int64}
        @test norm_v1 == [1,0,1,0,1,0]

        @test norm_v2 isa Vector{Int64}
        @test norm_v2 == [1,0,0,1,0,0]

        @test gen_v isa Vector{Int64}
        @test gen_v == [1,0,0,1,0,0,1,0,0]

        @test_throws DomainError convert(Vector{Int64}, s1, rep="no-signaling")
    end

    @testset "Vector{Int64} -> DeterministicStrategy" begin
        norm_v1 = [1,0,0,1,0,0,1,0,0]
        gen_v1 = [1,0,0,0,1,0,0,0,1,0,0,0]

        norm_v2 = [0,0,1,0,0,0,0,0,0]

        scenario = LocalSignaling(3,4,2)

        norm_s1 = convert(DeterministicStrategy, norm_v1, scenario, rep="normalized")
        norm_s2 = convert(DeterministicStrategy, norm_v2, scenario, rep="normalized")

        gen_s1 = convert(DeterministicStrategy, gen_v1, scenario, rep="generalized")

        @test norm_s1 isa DeterministicStrategy
        @test norm_s2 isa DeterministicStrategy
        @test gen_s1 isa DeterministicStrategy

        @test norm_s1 == [1 1 1;0 0 0;0 0 0;0 0 0]
        @test norm_s2 == [0 0 0;0 0 0;1 0 0;0 1 1]
        @test gen_s1 == [1 1 1;0 0 0;0 0 0;0 0 0]

        @test_throws DomainError convert(DeterministicStrategy, norm_v1, scenario, rep="no-signaling")
    end
end

@testset "deterministic_strategies(scenario::BlackBox)" begin
    @test deterministic_strategies(BlackBox(1,5)) == [[1 1 1 1 1]]
    @test deterministic_strategies(BlackBox(2,2)) == [
        [1 1;0 0],[0 1;1 0],[1 0;0 1],[0 0;1 1]
    ]
    @test deterministic_strategies(BlackBox(2,3)) == [
        [1 1 1;0 0 0],[0 1 1;1 0 0],[1 0 1;0 1 0],[0 0 1;1 1 0],
        [1 1 0;0 0 1],[0 1 0;1 0 1],[1 0 0;0 1 1],[0 0 0;1 1 1]
    ]

    @testset "spot checks" begin
        strategies = deterministic_strategies(BlackBox(5,7))

        @test length(strategies) == 5^7
        @test length(unique(strategies)) == 5^7
        @test all(s -> is_deterministic(s), strategies)
    end
end

end
