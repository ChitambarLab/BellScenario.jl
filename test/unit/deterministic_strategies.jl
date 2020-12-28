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
