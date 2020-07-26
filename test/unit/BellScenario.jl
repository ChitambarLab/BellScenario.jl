using Test

@testset "./src/BellScenario.jl" begin

using BellScenario
using QBase: QMath

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

    @testset "Strategy multiplication" begin
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

    @testset "Strategy DomainError" begin
        scenario_error = BlackBox(4,3)

        @test_throws DomainError Strategy([1 0;0 1], scenario_error)
        @test_throws DomainError Strategy([0.5 0.5;1 0.5])
    end
end

@testset "strategy_dims()" begin
    @testset "black box scenario" begin
        @test strategy_dims(BlackBox(2,5)) == (5,2)
    end

    @testset "bipartite scenario" begin
        @test strategy_dims(Bipartite((3,1),(1,4))) == (4,3)
        @test strategy_dims(Bipartite((2,2),(2,2))) == (4,4)
    end
end


@testset "Game()" begin
    @testset "int valued game" begin
        G = Game([1 0 0;0 1 0;0 0 1], 2)

        @test G isa Game{Int64}
        @test G == [1 0 0;0 1 0;0 0 1]
        @test G.game == [1 0 0;0 1 0;0 0 1]
        @test G.β == 2
    end

    @testset "float valued game" begin
        G = Game([0.5 0 0;0 0.5 0;0 0 0.5], 1.0)

        @test G isa Game{Float64}
        @test G == 1/2*[1 0 0;0 1 0;0 0 1]
        @test G.game == 1/2*[1 0 0;0 1 0;0 0 1]
        @test G.β == 1
    end

    @testset "errors" begin
        @test_throws MethodError Game([1 0 0;0 1 0;0 0 1], 1im)
    end
end

@testset "BellGame()" begin
    BG = BellGame([1 0 0;0 1 0;0 0 1], 2)

    @test BG.β == 2
    @test BG == [1 0 0;0 1 0;0 0 1]

    @test BG isa AbstractGame{Int64}
    @test BG isa BellGame
end

end
