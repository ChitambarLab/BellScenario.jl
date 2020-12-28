using Test

@testset "./src/scenarios.jl" begin

using BellScenario

@testset "BlackBox()" begin
    BB = BlackBox(2,3)

    @test BB.num_in == 3
    @test BB.num_out == 2

    @test BB isa BlackBox

    @test_throws DomainError BlackBox(0,2)
    @test_throws DomainError BlackBox(2,0)
    @test_throws DomainError BlackBox(0,-1)
end

@testset "BipartiteNonSignaling()" begin
    scenario = BipartiteNonSignaling(1,2,3,4)
    @test scenario isa Scenario
    @test scenario isa BipartiteNonSignaling

    @test scenario.A == 1
    @test scenario.B == 2
    @test scenario.X == 3
    @test scenario.Y == 4

    @testset "DomainErrors" begin
        @test_throws DomainError BipartiteNonSignaling(0,2,2,2)
        @test_throws DomainError BipartiteNonSignaling(2,0,2,2)
        @test_throws DomainError BipartiteNonSignaling(2,2,0,2)
        @test_throws DomainError BipartiteNonSignaling(2,2,2,0)
    end
end

@testset "Bipartite()" begin
    @testset "calling w/ BlackBox parameters, no Comm" begin
        BS = Bipartite((2,3),(4,5))

        @test BS isa Bipartite
        @test BS isa Scenario

        @test BS.A == BlackBox(2,3)
        @test BS.B == BlackBox(4,5)
        @test BS.dits == 1
        @test BS.bidirectional == false
    end

    @testset "calling w/ BlackBox parameters, no Comm" begin
        BS = Bipartite((2,2),(2,2), dits =3, bidirectional=true)

        @test BS isa Bipartite
        @test BS isa Scenario

        @test BS.A == BlackBox(2,2)
        @test BS.B == BlackBox(2,2)
        @test BS.dits == 3
        @test BS.bidirectional == true
    end

    @test_throws DomainError Bipartite((2,3),(4,5),dits=0)
end

@testset "LocalSignaling()" begin
    @testset "construction" begin
        scenario = LocalSignaling(3,4,2)

        @test scenario isa Scenario
        @test scenario isa LocalSignaling


        @test scenario.X == 3
        @test scenario.Y == 4
        @test scenario.d ==  2
    end

    @testset "DomainError" begin
        @test_throws DomainError LocalSignaling(0,2,2)
    end
end
end
