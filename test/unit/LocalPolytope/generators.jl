using Test

@testset "./sr/LocalPolytope/generaotrs.jl" begin

using BellScenario

@testset "generator_vertex()" begin
    @testset "simple 3-2-3 example" begin
        PM = PrepareAndMeasure(3,3,2)

        D = DeterministicStrategy([0 0 0;1 0 1;0 1 0],PM)

        generator = LocalPolytope.generator_vertex(D, PM)

        @test generator == [1 1 0;0 0 1;0 0 0]
        @test generator isa DeterministicStrategy
    end

    @testset "simple 6-4-4 example" begin
        PM = PrepareAndMeasure(6,4,4)
        D = DeterministicStrategy([1 0 1 1 0 0;0 0 0 0 1 0;0 1 0 0 0 0;0 0 0 0 0 1], PM)

        generator = LocalPolytope.generator_vertex(D, PM)

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
        PM = PrepareAndMeasure(6,6,2)

        @test LocalPolytope.generator_facet(BG, PM) == [1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1]
    end

    @testset "ambiguous guessing games" begin
        BG = BellGame([2 0 0;0 2 0;0 0 2;1 1 1],4)
        PM = PrepareAndMeasure(3,4,2)

        @test LocalPolytope.generator_facet(BG, PM) == [2 0 0;1 1 1;0 2 0;0 0 2]

        BG = BellGame([1 1 1 1 0;0 0 3 0 0;0 0 0 3 0;0 1 0 0 2;3 0 0 0 0], 6)
        PM = PrepareAndMeasure(5,5,2)

        @test LocalPolytope.generator_facet(BG,  PM) == [3 0 0 0 0;1 1 1 1 0;0 3 0 0 0;0 0 3 0 0;0 0 0 1 2]
    end

    @testset "error guessing games" begin
        BG = BellGame([1 0 0 0;0 1 0 1;0 0 1 1;0 1 1 0], 3)
        PM = PrepareAndMeasure(4,4,2)

        @test LocalPolytope.generator_facet(BG, PM) == [1 1 0 0;1 0 1 0;0 1 1 0;0 0 0 1]
    end

    @testset "coarse-grained output" begin
        BG = BellGame([0 1 0;1 0 0;0 0 1;0 1 0],2)
        PM = PrepareAndMeasure(3,4,2)

        @time @test LocalPolytope.generator_facet(BG, PM) == [1 0 0;1 0 0;0 1 0;0 0 1]
    end

    @testset "7x7 case raw permutations fails here" begin
        BG = BellGame([
                1 1 1 1 1 1 0 0;0 0 0 0 0 0 0 1;
                0 1 1 1 1 1 1 0;1 0 1 1 1 1 1 0;
                1 1 1 1 0 1 1 0;1 1 0 1 1 1 1 0;
                1 1 1 0 1 1 1 0;1 1 1 1 1 0 1 0;
            ], 8)

        PM = PrepareAndMeasure(8,8,2)

        @test LocalPolytope.generator_facet(BG, PM) == [
            1 1 1 1 1 1 0 0;1 1 1 1 1 0 1 0;
            1 1 1 1 0 1 1 0;1 1 1 0 1 1 1 0;
            1 1 0 1 1 1 1 0;1 0 1 1 1 1 1 0;
            0 1 1 1 1 1 1 0;0 0 0 0 0 0 0 1;
        ]
    end

    @testset "testing cases with duplicate columns" begin
        PM = PrepareAndMeasure(5,5,2)
        BG1 = BellGame([1 0 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 1 0;0 0 0 0 1], 2)
        BG2 = BellGame([1 0 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1;0 0 0 0 1], 2)
        BG3 = BellGame([0 0 0 0 1;0 0 1 0 0;0 0 1 0 0;0 0 0 1 0;0 1 0 0 0], 2)

        match = [1 0 0 0 0;1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0]

        @test LocalPolytope.generator_facet(BG1, PM) == match
        @test LocalPolytope.generator_facet(BG2, PM) == match
        @test LocalPolytope.generator_facet(BG3, PM) == match
    end

    @testset "test case where values are re-mapped" begin
        BG = BellGame([4 2 2 0 0;
                4 0 0 2 2; 3 0 3 3 0;
                0 3 0 3 0; 0 0 4 0 2
            ], 12)
        PM = PrepareAndMeasure(5,5,2)

        @test LocalPolytope.generator_facet(BG, PM) == [4 2 2 0 0;4 0 0 2 2;3 3 0 3 0;0 3 0 0 3;0 0 2 4 0]
    end

    @testset "historic failure no longer fails" begin
        BG = BellGame([1 1 0 0 0; 1 0 1 0 0; 0 1 1 0 0; 0 0 0 0 1; 0 0 0 0 1], 3)
        PM = PrepareAndMeasure(5,5,2)

        @test LocalPolytope.generator_facet(BG,PM)  == [1 1 0 0 0; 1 0 1 0 0; 0 1 1 0 0; 0 0 0 1 0; 0 0 0 1 0]
    end
end

end
