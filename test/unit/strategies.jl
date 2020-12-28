using Test

@testset "./src/strategies.jl" begin

using BellScenario

@testset "Strategy()" begin
    @testset "with scenario" begin
    probs = [1 0 0;0 1 0;0 0 1]

        scenario = BlackBox(3,3)

        S = Strategy(probs, scenario)

        @test S == probs
        @test S isa Strategy
        @test S.scenario isa BlackBox
        @test S.scenario == scenario
    end

    @testset "without scenario, initializes as BlackBox" begin
        S = Strategy([0.5 0.5 0.5;0.5 0.5 0.5])

        @test S == 1/2*[1 1 1;1 1 1]

        @test S.scenario == BlackBox(2,3)
    end

    @testset "Strategy DomainError" begin
        scenario_error = BlackBox(3,4)

        @test_throws DomainError Strategy([1 0;0 1], scenario_error)
        @test_throws DomainError Strategy([0.5 0.5;1 0.5])
    end
end

@testset "Strategy multiplication" begin
    @testset "Strategy * Strategy" begin
        S1 = Strategy([1 0;0 1;0 0])
        S2 = Strategy([1 0 1;0 1 0])

        S_prod1 = S1*S2
        S_prod2 = S2*S1

        @test S_prod1 == [1 0 1;0 1 0;0 0 0]
        @test S_prod1 isa Strategy
        @test S_prod1.scenario == BlackBox(3,3)

        @test S_prod2 == [1 0;0 1]
        @test S_prod2 isa Strategy
        @test S_prod2.scenario == BlackBox(2,2)

        scenario = BipartiteSignaling((3,1),(1,3))
        S_prod = *(S1, S2, scenario)

        @test S_prod == [1 0 1;0 1 0;0 0 0]
        @test S_prod isa Strategy
        @test S_prod.scenario == BipartiteSignaling((3,1),(1,3))

        @test_throws DomainError (*(S2, S1, scenario))
    end

    @testset "DeterministicStrategy * DeterministicStrategy" begin
        S1 = DeterministicStrategy([1 0;0 1;0 0])
        S2 = DeterministicStrategy([1 0 1;0 1 0])

        S_prod1 = S1*S2
        S_prod2 = S2*S1

        @test S_prod1 == [1 0 1;0 1 0;0 0 0]
        @test S_prod1 isa DeterministicStrategy
        @test S_prod1.scenario == BlackBox(3,3)

        @test S_prod2 == [1 0;0 1]
        @test S_prod2 isa DeterministicStrategy
        @test S_prod2.scenario == BlackBox(2,2)

        scenario = BipartiteSignaling((3,1),(1,3))
        S_prod = *(S1, S2, scenario)

        @test S_prod == [1 0 1;0 1 0;0 0 0]
        @test S_prod isa DeterministicStrategy
        @test S_prod.scenario == BipartiteSignaling((3,1),(1,3))

        @test_throws DomainError (*(S2, S1, scenario))
    end

    @testset "DeterministicStrategy * Strategy is stochastic" begin
        S1 = Strategy([1 0;0 1;0 0])
        S2 = DeterministicStrategy([1 0 1;0 1 0])

        S_prod1 = S1*S2
        S_prod2 = S2*S1

        @test S_prod1 == [1 0 1;0 1 0;0 0 0]
        @test S_prod1 isa Strategy
        @test S_prod1.scenario == BlackBox(3,3)

        @test S_prod2 == [1 0;0 1]
        @test S_prod2 isa Strategy
        @test S_prod2.scenario == BlackBox(2,2)

        scenario = BipartiteSignaling((3,1),(1,3))
        S_prod = *(S1, S2, scenario)

        @test S_prod == [1 0 1;0 1 0;0 0 0]
        @test S_prod isa Strategy
        @test S_prod.scenario == BipartiteSignaling((3,1),(1,3))

        @test_throws DomainError (*(S2, S1, scenario))
    end
end

@testset "random_strategy()" begin
    @testset "simple cases"  begin
        strat = random_strategy(3,3)

        @test strat isa Strategy
        @test size(strat) == (3,3)

        strat = random_strategy(9,7)

        @test strat isa Strategy
        @test size(strat) == (7,9)
    end
end

@testset "strategy_dims()" begin
    @testset "black box scenario" begin
        @test strategy_dims(BlackBox(5,2)) == (5,2)
    end

    @testset "prepare and measure" begin
        @test strategy_dims(LocalSignaling(3,5,2)) == (5,3)
    end

    @testset "bipartite non-signaling scenario" begin
        @test strategy_dims(BipartiteNonSignaling(2,2,2,2)) == (4,4)
        @test strategy_dims(BipartiteNonSignaling(2,3,4,5)) == (6,20)
    end

    @testset "bipartite scenario" begin
        @test strategy_dims(BipartiteSignaling((3,1),(1,4))) == (3,4)
        @test strategy_dims(BipartiteSignaling((2,2),(2,2))) == (4,4)
    end
end

end
