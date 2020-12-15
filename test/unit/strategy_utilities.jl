using Test

@testset "./src/strategy_utilities.jl" begin

using BellScenario

@testset "convert(DeterministicStrategy, ...)" begin
    @testset "chsh vertices no-signaling" begin
        scenario = BipartiteNoSignaling(2,2,2,2)
        vertex = [1,0,1,0,1,0,0,0]

        det_strat = convert(DeterministicStrategy, vertex, scenario)

        @test det_strat isa DeterministicStrategy
        @test det_strat == [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]

        strat = convert(Strategy, vertex, scenario)

        @test strat isa Strategy
        @test strat == [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]
    end

    @testset "chsh vertices normalized" begin
        scenario = BipartiteNoSignaling(2,2,2,2)
        vertex = [1,0,0,0,1,0,1,0,0,0,1,0]

        det_strat = convert(DeterministicStrategy, vertex, scenario, rep="normalized")

        @test det_strat isa DeterministicStrategy
        @test det_strat == [1 0 1 0;0 1 0 1;0 0 0 0;0 0 0 0]

        strat = convert(Strategy, vertex, scenario, rep="normalized")

        @test strat isa Strategy
        @test strat == [1 0 1 0;0 1 0 1;0 0 0 0;0 0 0 0]
    end

    @testset "chsh vertices generalized" begin
        scenario = BipartiteNoSignaling(2,2,2,2)
        vertex = [1,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0]

        det_strat = convert(DeterministicStrategy, vertex, scenario, rep="generalized")

        @test det_strat isa DeterministicStrategy
        @test det_strat == [1 0 1 0;0 1 0 1;0 0 0 0;0 0 0 0]

        strat = convert(Strategy, vertex, scenario, rep="generalized")

        @test strat isa Strategy
        @test strat == [1 0 1 0;0 1 0 1;0 0 0 0;0 0 0 0]
    end

    @testset "non-signaling vertices" begin
        scenario = BipartiteNoSignaling(2,2,2,2)

        ns_vertex1 = [1/2,1/2,1/2,1/2,0,0,0,1/2]
        ns_strat1 = convert(Strategy, ns_vertex1, scenario)

        @test ns_strat1 isa Strategy
        @test ns_strat1 == 0.5*[0 0 0 1;1 1 1 0;1 1 1 0;0 0 0 1]

        ns_vertex2 = [1/2,1/2,1/2,1/2,1/2,0,1/2,1/2]
        ns_strat2 = convert(Strategy, ns_vertex2, scenario)

        @test ns_strat2 isa Strategy
        @test ns_strat2 == 0.5*[1 0 1 1;0 1 0 0;0 1 0 0;1 0 1 1]
    end

    @testset "DomainError strategy not deterministic" begin
        scenario = BipartiteNoSignaling(2,2,2,2)
        behavior = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]

        @test_throws DomainError convert(DeterministicStrategy, behavior, scenario)
    end

    @testset "DomainError not valid rep" begin
        scenario = BipartiteNoSignaling(2,2,2,2)
        behavior = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
        
        @test_throws DomainError convert(Strategy, behavior, scenario, rep="fixed-direction")
    end
end

end
