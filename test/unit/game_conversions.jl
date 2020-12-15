using Test
using XPORTA: IEQ

@testset "./src/game_conversions.jl" begin

using BellScenario

@testset "Conversions to BellGames" begin
    @testset  "Vector{Int64} -> BellGame, LocalSignaling" begin
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

    @testset "IEQ -> BellGames, BlackBox" begin
        gen_ieq = IEQ(inequalities = [1 1 0 0 0 1 2;1 1 0 0 0 0 1] )
        scenario = BlackBox(3,2)

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

    @testset "Vector{Int64} -> BellGame, BipartiteNoSignaling" begin
        scenario = BipartiteNoSignaling(2,2,2,2)

        ns_facet = [1, 0, 1, 0, -1, -1, -1, 1, 1]
        bell_game = convert(BellGame, ns_facet, scenario)

        @test bell_game isa BellGame
        @test bell_game == [1 0 0 1;1 1 1 0;1 1 1 0;0 1 1 0]
        @test bell_game.β == 3

        norm_facet = [0, 0, 1, 0, 1, -1, 0, 0, 0, 1, 0, 1, 1]
        bell_game = convert(BellGame, norm_facet, scenario, rep="normalized")

        @test bell_game isa BellGame
        @test bell_game.β == 2
        @test bell_game == [0 1 0 1;0 2 0 0;1 0 0 1;0 1 0 0]
    end
end

@testset "Conversions from BellGames" begin

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

    @testset "BellGame -> Facet BipartiteNoSignaling" begin
        @testset "chsh spot check" begin
            BG = BellGame([1 0 0 1;1 1 1 0;1 1 1 0;0 1 1 0], 3)
            scenario = BipartiteNoSignaling(2,2,2,2)

            facet = convert(Vector{Int64}, BG, scenario)

            @test facet isa Vector{Int64}
            @test length(facet) == LocalPolytope.vertex_dims(scenario, "no-signaling") + 1
            @test facet == [1,0,1,0,-1,-1,-1,1,1]
        end
    end
end

@testset "identity conversions" begin
    @testset "BellGames -> IEQ -> BellGames" begin
        bell_games = BellGame.([[2 0 0 0;0 2 0 0;0 0 2 0;1 1 1 0],[1 0 0 0;1 0 0 0;1 0 0 0;0 0 0 0]],[4,1])
        scenario = LocalSignaling(4,4,2)

        ieq = convert(IEQ, bell_games)
        converted_bell_games = convert(Vector{BellGame}, ieq, scenario)

        @test converted_bell_games == bell_games
    end

    @testset "BipartiteNoSignaling facet -> BellGame -> facet" begin
        @testset "simple scenarios full check" for scenario_params in [
            (2,2,2,2),(3,2,2,2),(2,3,2,2),(2,2,3,2),(2,2,2,3),
            (3,3,2,2),(3,2,3,2),(3,2,2,3),(2,3,3,2),(2,3,2,3),
            (2,2,3,3)
        ]
            scenario = BipartiteNoSignaling(scenario_params...)

            vertices = LocalPolytope.vertices(scenario)
            facet_dict = LocalPolytope.facets(vertices)

            for facet in facet_dict["facets"]
                BG = convert(BellGame, facet, scenario)

                @test BG isa BellGame
                println(BG)

                facet2 = convert(Vector{Int64}, BG, scenario)

                @test facet2 isa Vector{Int64}

                @test facet2 == facet
            end
        end
    end
end

end
