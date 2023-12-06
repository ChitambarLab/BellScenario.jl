using Test, LinearAlgebra, Polyhedra, XPORTA

@testset "test/unit/LocalPolytope.jl" begin

using BellScenario

@testset "running all tests for module LocalPolytope.jl" begin

    @testset "LocalPolytope.vrep()" begin
        @testset "no arguments" begin
            scenario = LocalSignaling(3,3,2)
            local_poly = LocalPolytope.vrep(scenario)

            @test local_poly isa XPORTA.Polyhedron
            @test npoints(local_poly) == 21

            vertices = map(p -> convert.(Int64,p), points(local_poly))
            @test LocalPolytope.vertices(scenario) == vertices
        end

        @testset "passing arguments" begin
            scenario = LocalSignaling(3,3,2)
            local_poly = LocalPolytope.vrep(scenario, rank_d_only=true)

            @test local_poly isa XPORTA.Polyhedron
            @test npoints(local_poly) == 18

            vertices = map(p -> convert.(Int64,p), points(local_poly))
            @test LocalPolytope.vertices(scenario, rank_d_only=true) == vertices
        end
    end

    println("running LocalPolytope.jl unit tests:")
    for test in readdir("./test/unit/LocalPolytope/")
    # run only julia files in test directory
    if occursin(r"^.*\.jl$", test)
        println("./unit/LocalPolytope/$test")
        @time include("./LocalPolytope/$test")
    end
    end

end

end
