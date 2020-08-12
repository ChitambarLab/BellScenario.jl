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

        @test S.scenario == BlackBox(3,2)
    end

    @testset "Strategy DomainError" begin
        scenario_error = BlackBox(4,3)

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

        scenario = Bipartite((3,1),(1,3))
        S_prod = *(S1, S2, scenario)

        @test S_prod == [1 0 1;0 1 0;0 0 0]
        @test S_prod isa Strategy
        @test S_prod.scenario == Bipartite((3,1),(1,3))

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

        scenario = Bipartite((3,1),(1,3))
        S_prod = *(S1, S2, scenario)

        @test S_prod == [1 0 1;0 1 0;0 0 0]
        @test S_prod isa DeterministicStrategy
        @test S_prod.scenario == Bipartite((3,1),(1,3))

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

        scenario = Bipartite((3,1),(1,3))
        S_prod = *(S1, S2, scenario)

        @test S_prod == [1 0 1;0 1 0;0 0 0]
        @test S_prod isa Strategy
        @test S_prod.scenario == Bipartite((3,1),(1,3))

        @test_throws DomainError (*(S2, S1, scenario))
    end
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

        PM = PrepareAndMeasure(3,4,2)

        norm_s1 = convert(DeterministicStrategy, norm_v1, PM, rep="normalized")
        norm_s2 = convert(DeterministicStrategy, norm_v2, PM, rep="normalized")

        gen_s1 = convert(DeterministicStrategy, gen_v1, PM, rep="generalized")

        @test norm_s1 isa DeterministicStrategy
        @test norm_s2 isa DeterministicStrategy
        @test gen_s1 isa DeterministicStrategy

        @test norm_s1 == [1 1 1;0 0 0;0 0 0;0 0 0]
        @test norm_s2 == [0 0 0;0 0 0;1 0 0;0 1 1]
        @test gen_s1 == [1 1 1;0 0 0;0 0 0;0 0 0]

        @test_throws DomainError convert(DeterministicStrategy, norm_v1, PM, rep="no-signaling")
    end
end

@testset "strategy_dims()" begin
    @testset "black box scenario" begin
        @test strategy_dims(BlackBox(2,5)) == (5,2)
    end

    @testset "prepare and measure" begin
        @test strategy_dims(PrepareAndMeasure(3,5,2)) == (5,3)
    end

    @testset "bipartite scenario" begin
        @test strategy_dims(Bipartite((3,1),(1,4))) == (4,3)
        @test strategy_dims(Bipartite((2,2),(2,2))) == (4,4)
    end
end

@testset "DeterministicStrategy()" begin
    D = DeterministicStrategy([1.0 0.0;0.0 1.0;0 0])

    @test D isa DeterministicStrategy
    @test D.conditionals isa Matrix{Int64}
    @test D == [1 0;0 1;0 0]
    @test D.scenario == BlackBox(2,3)

    @test_throws DomainError DeterministicStrategy([0.5 0.5;0.5 0.5])
    @test_throws DomainError DeterministicStrategy([1 0;0 1], BlackBox(3,2))
end

@testset "is_deterministic()" begin
    @test is_deterministic([1 0;0 1])
    @test !is_deterministic([0.5 0.5;0.5 0.5])
    @test is_deterministic([1.0 1.0;0.0 0.0])
end

end
