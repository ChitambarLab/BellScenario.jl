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

@testset "BellGame conversions" begin
    @testset "BellGame -> Vector{Int64}" begin
        BG = BellGame([2 0 0;1 1 1;0 2 0;0 0 2], 4)

        gen_facet = convert(Vector{Int64}, BG, rep="generalized")

        @test gen_facet isa Vector{Int64}
        @test gen_facet == [2,1,0,0,0,1,2,0,0,1,0,2,4]

        norm_facet = convert(Vector{Int64}, BG, rep="normalized")

        @test norm_facet isa Vector{Int64}
        @test norm_facet == [2,1,0,0,1,2,-2,-1,-2,2]

        @test_throws DomainError convert(Vector{Int64}, BG, rep="no-signaling")
    end

    @testset  "Vector{Int64} -> BellGame" begin
        norm_facet = [2,1,0,0,1,2,-2,-1,-2,2]
        gen_facet = [2,1,0,0,0,1,2,0,0,1,0,2,4]
        scenario = LocalSignaling(3,4,2)

        bell_game = convert(BellGame, norm_facet, scenario, rep="normalized")
        bell_game2 = convert(BellGame, gen_facet, scenario, rep="generalized")

        @test bell_game isa BellGame
        @test bell_game2 isa BellGame
        @test bell_game == [2 0 0;1 1 1;0 2 0;0 0 2]
        @test bell_game.β == 4
        @test bell_game2 == [2 0 0;1 1 1;0 2 0;0 0 2]
        @test bell_game2.β == 4

        @test_throws DomainError convert(BellGame, norm_facet, scenario, rep="no-signaling")
    end

    @testset "BellGames -> IEQ" begin
        bell_games = BellGame.([[1 0;0 1],[1 0;0 0]],[2,1])
        gen_ieq = convert(IEQ, bell_games, rep="generalized")

        @test gen_ieq isa IEQ
        @test gen_ieq.inequalities == [1 0 0 1 2;1 0 0 0 1]

        norm_ieq = convert(IEQ, bell_games, rep="normalized")

        @test norm_ieq isa IEQ
        @test norm_ieq.inequalities == [1 -1 1;1 0 1]

        norm_ieq2 = convert(IEQ, bell_games)

        @test norm_ieq2 isa IEQ
        @test norm_ieq2.inequalities == [1 -1 1;1 0 1]
    end

    @testset "IEQ -> BellGames" begin
        gen_ieq = IEQ(inequalities = [1 1 0 0 0 1 2;1 1 0 0 0 0 1] )
        scenario = BlackBox(2,3)

        bell_games = convert(Vector{BellGame}, gen_ieq, scenario, rep="generalized")

        @test bell_games isa Vector{BellGame}
        @test bell_games == BellGame.([[1 0;1 0;0 1],[1 0;1 0;0 0]],[2,1])

        norm_ieq = IEQ(inequalities = [1 1 -1 -1 1;1 1 0 0 1])

        bell_games = convert(Vector{BellGame}, norm_ieq, scenario, rep="normalized")
        bell_games2 = convert(Vector{BellGame}, norm_ieq, scenario)

        @test bell_games isa Vector{BellGame}
        @test bell_games2 isa Vector{BellGame}
        @test bell_games == BellGame.([[1 0;1 0;0 1],[1 0;1 0;0 0]],[2,1])
        @test bell_games2 == BellGame.([[1 0;1 0;0 1],[1 0;1 0;0 0]],[2,1])
    end

    @testset "BellGames -> IEQ -> BellGames" begin
        bell_games = BellGame.([[2 0 0 0;0 2 0 0;0 0 2 0;1 1 1 0],[1 0 0 0;1 0 0 0;1 0 0 0;0 0 0 0]],[4,1])
        scenario = LocalSignaling(4,4,2)

        ieq = convert(IEQ, bell_games)
        converted_bell_games = convert(Vector{BellGame}, ieq, scenario)

        @test converted_bell_games == bell_games
    end
end

end
