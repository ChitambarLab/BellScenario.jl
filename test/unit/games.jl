using Test

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

end
