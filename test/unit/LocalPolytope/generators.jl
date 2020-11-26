using Test

@testset "./sr/LocalPolytope/generaotrs.jl" begin

using BellScenario

@testset "generator_vertex()" begin
    @testset "simple 3-2-3 example" begin
        scenario = LocalSignaling(3,3,2)

        D = DeterministicStrategy([0 0 0;1 0 1;0 1 0],scenario)

        generator = LocalPolytope.generator_vertex(D, scenario)

        @test generator == [1 1 0;0 0 1;0 0 0]
        @test generator isa DeterministicStrategy
    end

    @testset "simple 6-4-4 example" begin
        scenario = LocalSignaling(6,4,4)
        D = DeterministicStrategy([1 0 1 1 0 0;0 0 0 0 1 0;0 1 0 0 0 0;0 0 0 0 0 1], scenario)

        generator = LocalPolytope.generator_vertex(D, scenario)

        @test generator == [1 1 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1]
        @test D isa DeterministicStrategy
    end
end

@testset "_sort_perm_groups()" begin
    @test LocalPolytope._sort_perm_groups(
            [(1,1),(3,2),(4,3),(2,4),(5,5)],
            [[1],[2,3],[4,5]]
        ) == [(1,1),(4,3),(3,2),(5,5),(2,4)]
end

@testset "_perm_increase_lexico_score()" begin
    game = [1 1 1 0;0 0 2 0;0 2 0 0;2 0 0 0]
    allowed_col_perms = [collect(1:4)]
    target_row = 1

    canonical_game = LocalPolytope._perm_increase_lexico_score(game, target_row, allowed_col_perms)

    @test canonical_game == [2 0 0 0;1 1 1 0;0 2 0 0;0 0 2 0]
end

@testset "generator_facet()" begin
    @testset "success guessing game" begin
        BG = BellGame([0 1 0 0 0 0;0 0 0 0 1 0;1 0 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 0 1],2)
        scenario = LocalSignaling(6,6,2)

        @test LocalPolytope.generator_facet(BG, scenario) == [1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1]
    end

    @testset "ambiguous guessing games" begin
        BG = BellGame([2 0 0;0 2 0;0 0 2;1 1 1],4)
        scenario = LocalSignaling(3,4,2)

        @test LocalPolytope.generator_facet(BG, scenario) == [2 0 0;1 1 1;0 2 0;0 0 2]

        BG = BellGame([1 1 1 1 0;0 0 3 0 0;0 0 0 3 0;0 1 0 0 2;3 0 0 0 0], 6)
        scenario = LocalSignaling(5,5,2)

        @test LocalPolytope.generator_facet(BG,  scenario) == [3 0 0 0 0;1 1 1 1 0;0 3 0 0 0;0 0 3 0 0;0 0 0 1 2]
    end

    @testset "error guessing games" begin
        BG = BellGame([1 0 0 0;0 1 0 1;0 0 1 1;0 1 1 0], 3)
        scenario = LocalSignaling(4,4,2)

        @test LocalPolytope.generator_facet(BG, scenario) == [1 1 0 0;1 0 1 0;0 1 1 0;0 0 0 1]
    end

    @testset "coarse-grained output" begin
        BG = BellGame([0 1 0;1 0 0;0 0 1;0 1 0],2)
        scenario = LocalSignaling(3,4,2)

        @time @test LocalPolytope.generator_facet(BG, scenario) == [1 0 0;1 0 0;0 1 0;0 0 1]
    end

    @testset "7x7 case raw permutations fails here" begin
        BG = BellGame([
                1 1 1 1 1 1 0 0;0 0 0 0 0 0 0 1;
                0 1 1 1 1 1 1 0;1 0 1 1 1 1 1 0;
                1 1 1 1 0 1 1 0;1 1 0 1 1 1 1 0;
                1 1 1 0 1 1 1 0;1 1 1 1 1 0 1 0;
            ], 8)

        scenario = LocalSignaling(8,8,2)

        @test LocalPolytope.generator_facet(BG, scenario) == [
            1 1 1 1 1 1 0 0;1 1 1 1 1 0 1 0;
            1 1 1 1 0 1 1 0;1 1 1 0 1 1 1 0;
            1 1 0 1 1 1 1 0;1 0 1 1 1 1 1 0;
            0 1 1 1 1 1 1 0;0 0 0 0 0 0 0 1;
        ]
    end

    @testset "testing cases with duplicate columns" begin
        scenario = LocalSignaling(5,5,2)
        BG1 = BellGame([1 0 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 1 0;0 0 0 0 1], 2)
        BG2 = BellGame([1 0 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1;0 0 0 0 1], 2)
        BG3 = BellGame([0 0 0 0 1;0 0 1 0 0;0 0 1 0 0;0 0 0 1 0;0 1 0 0 0], 2)

        match = [1 0 0 0 0;1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0]

        @test LocalPolytope.generator_facet(BG1, scenario) == match
        @test LocalPolytope.generator_facet(BG2, scenario) == match
        @test LocalPolytope.generator_facet(BG3, scenario) == match
    end

    @testset "test case where values are re-mapped" begin
        BG = BellGame([4 2 2 0 0;
                4 0 0 2 2; 3 0 3 3 0;
                0 3 0 3 0; 0 0 4 0 2
            ], 12)
        scenario = LocalSignaling(5,5,2)

        @test LocalPolytope.generator_facet(BG, scenario) == [4 2 2 0 0;4 0 0 2 2;3 3 0 3 0;0 3 0 0 3;0 0 2 4 0]
    end

    @testset "historic failure no longer fails" begin
        BG = BellGame([1 1 0 0 0; 1 0 1 0 0; 0 1 1 0 0; 0 0 0 0 1; 0 0 0 0 1], 3)
        scenario = LocalSignaling(5,5,2)

        @test LocalPolytope.generator_facet(BG,scenario)  == [1 1 0 0 0; 1 0 1 0 0; 0 1 1 0 0; 0 0 0 1 0; 0 0 0 1 0]
    end
end

end
