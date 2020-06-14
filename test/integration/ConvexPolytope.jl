using Test

@testset "src/ConvexPolytope.jl" begin

using BellComm: ConvexPolytope

@testset "facet_constraints()" begin
    dir = "./test/integration/files/"

    @testset "cleanup = false" begin
        test_dir = mkpath(dir * "test_porta_tmp")

        @test_throws SystemError readdir(test_dir * "porta_tmp")

        poi_file = test_dir * "/facet_constraints_test.poi"
        ieq_file = test_dir * "/facet_constraints_test.poi.ieq"

        @test !(isfile(poi_file))
        @test !(isfile(ieq_file))

        vertices = [[0;0;0],[1;0;0],[0;1;0],[1;1;1]]
        constraints = ConvexPolytope.facet_constraints(vertices, filename="facet_constraints_test", dir=test_dir, cleanup=false)
        @test constraints == ([[0 0 0 -1],[0 -1 0 1],[0 0 -1 1],[-1 1 1 -1]],[],[])

        @test isfile(poi_file)
        @test isfile(ieq_file)

        @test length(readdir(test_dir)) == 2

        rm(test_dir, recursive=true)
        @test_throws SystemError readdir(test_dir)
    end

    @testset "cleanup = true" begin

        vertices = [[0;0;0],[1;0;0],[0;1;0],[1;1;1]]
        constraints = ConvexPolytope.facet_constraints(vertices, dir=dir)
        @test constraints == (
            [[0 0 0 -1],[0 -1 0 1],[0 0 -1 1],[-1 1 1 -1]],
            [],
            []
        )

        @test_throws SystemError readdir(dir * "porta_tmp")
    end

    @testset "default filename used" begin
        test_dir = mkpath(dir * "test_porta_tmp")

        @test_throws SystemError readdir(test_dir * "porta_tmp")

        poi_file = test_dir * "/facet_constraints_tmp.poi"
        ieq_file = test_dir * "/facet_constraints_tmp.poi.ieq"

        @test !(isfile(poi_file))
        @test !(isfile(ieq_file))

        vertices = [[0;0;0],[1;0;0],[0;1;0],[1;1;1]]
        constraints = ConvexPolytope.facet_constraints(vertices, dir=test_dir, cleanup=false)
        @test constraints == (
            [[0 0 0 -1],[0 -1 0 1],[0 0 -1 1],[-1 1 1 -1]],
            [],
            []
        )

        @test isfile(poi_file)
        @test isfile(ieq_file)
        @test length(readdir(test_dir)) == 2

        rm(test_dir, recursive=true)
        @test_throws SystemError readdir(test_dir)
    end

    @testset "all defaults used" begin

        vertices = [[0;0;0],[1;0;0],[0;1;0],[1;1;1]]
        constraints = ConvexPolytope.facet_constraints(vertices)
        @test constraints == (
            [[0 0 0 -1],[0 -1 0 1],[0 0 -1 1],[-1 1 1 -1]],
            [],
            []
        )

        @test_throws SystemError readdir("./porta_tmp")
    end
end

end
