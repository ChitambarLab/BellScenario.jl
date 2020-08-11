using Test
using XPORTA: IEQ

@testset "./src/games.jl" begin

using BellScenario

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

@testset "Base.convert()" begin
    @testset "BellGames -> IEQ" begin
        bell_games = BellGame.([[1 0;0 1],[1 0;0 0]],[2,1])
        ieq = convert(IEQ, bell_games)

        @test ieq isa IEQ
        @test ieq.inequalities == [1 0 0 1 2;1 0 0 0 1]
    end

    @testset "IEQ -> BellGames" begin
        ieq = IEQ(inequalities = [1 0 1 0 0 1 2;1 0 1 0 0 0 1] )
        scenario = BlackBox(2,3)

        bell_games = convert(Vector{BellGame}, ieq, scenario)

        @test bell_games isa Vector{BellGame}
        @test bell_games == BellGame.([[1 0;1 0;0 1],[1 0;1 0;0 0]],[2,1])
    end

    @testset "BellGames -> IEQ -> BellGames" begin
        bell_games = BellGame.([[2 0 0 0;0 2 0 0;0 0 2 0;1 1 1 0],[1 0 0 0;1 0 0 0;1 0 0 0;0 0 0 0]],[4,1])
        scenario = PrepareAndMeasure(4,4,2)

        ieq = convert(IEQ, bell_games)
        converted_bell_games = convert(Vector{BellGame}, ieq, scenario)

        @test converted_bell_games == bell_games
    end
end

end
