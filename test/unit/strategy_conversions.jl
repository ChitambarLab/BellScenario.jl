using Test

@testset "./src/strategy_conversions.jl" begin

using BellScenario

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

        @test_throws DomainError convert(Vector{Int64}, s1, rep="non-signaling")
    end

    @testset "LocalSignaling :: Vector{Int64} -> DeterministicStrategy" begin
        norm_v1 = [1,0,0,1,0,0,1,0,0]
        gen_v1 = [1,0,0,0,1,0,0,0,1,0,0,0]

        norm_v2 = [0,0,1,0,0,0,0,0,0]

        scenario = LocalSignaling(3,4,2)

        norm_s1 = convert(DeterministicStrategy, norm_v1, scenario, rep="normalized")
        norm_s2 = convert(DeterministicStrategy, norm_v2, scenario, rep="normalized")

        gen_s1 = convert(DeterministicStrategy, gen_v1, scenario, rep="generalized")

        @test norm_s1 isa DeterministicStrategy
        @test norm_s2 isa DeterministicStrategy
        @test gen_s1 isa DeterministicStrategy

        @test norm_s1 == [1 1 1;0 0 0;0 0 0;0 0 0]
        @test norm_s2 == [0 0 0;0 0 0;1 0 0;0 1 1]
        @test gen_s1 == [1 1 1;0 0 0;0 0 0;0 0 0]

        @test_throws DomainError convert(DeterministicStrategy, norm_v1, scenario, rep="non-signaling")
    end
end

@testset "convert(DeterministicStrategy, ...)" begin
    @testset "chsh vertices non-signaling" begin
        scenario = BipartiteNonSignaling(2,2,2,2)
        vertex = [1,0,1,0,1,0,0,0]

        det_strat = convert(DeterministicStrategy, vertex, scenario)

        @test det_strat isa DeterministicStrategy
        @test det_strat == [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]

        strat = convert(Strategy, vertex, scenario)

        @test strat isa Strategy
        @test strat == [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]
    end

    @testset "chsh vertices normalized" begin
        scenario = BipartiteNonSignaling(2,2,2,2)
        vertex = [1,0,0,0,1,0,1,0,0,0,1,0]

        det_strat = convert(DeterministicStrategy, vertex, scenario, rep="normalized")

        @test det_strat isa DeterministicStrategy
        @test det_strat == [1 0 1 0;0 1 0 1;0 0 0 0;0 0 0 0]

        strat = convert(Strategy, vertex, scenario, rep="normalized")

        @test strat isa Strategy
        @test strat == [1 0 1 0;0 1 0 1;0 0 0 0;0 0 0 0]
    end

    @testset "chsh vertices generalized" begin
        scenario = BipartiteNonSignaling(2,2,2,2)
        vertex = [1,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0]

        det_strat = convert(DeterministicStrategy, vertex, scenario, rep="generalized")

        @test det_strat isa DeterministicStrategy
        @test det_strat == [1 0 1 0;0 1 0 1;0 0 0 0;0 0 0 0]

        strat = convert(Strategy, vertex, scenario, rep="generalized")

        @test strat isa Strategy
        @test strat == [1 0 1 0;0 1 0 1;0 0 0 0;0 0 0 0]
    end

    @testset "non-signaling vertices" begin
        scenario = BipartiteNonSignaling(2,2,2,2)

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
        scenario = BipartiteNonSignaling(2,2,2,2)
        behavior = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]

        @test_throws DomainError convert(DeterministicStrategy, behavior, scenario)
    end

    @testset "DomainError not valid rep" begin
        scenario = BipartiteNonSignaling(2,2,2,2)
        behavior = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]

        @test_throws DomainError convert(Strategy, behavior, scenario, rep="fixed-direction")
    end
end

end
