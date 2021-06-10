using Test, JSON

@testset "./src/LocalPolytope/adjacency_decomposition.jl" begin

using BellScenario

test_dir = "./test/integration/files/"

@testset "adjacent_facets" begin
    @testset "41-2-14 polytope" begin
        scenario = LocalSignaling(4,4,2)
        vertices = LocalPolytope.vertices(scenario)

        F = [1,0,0,0,1,0,0,0,1,-1,-1,-1,1]
        adj_facets = LocalPolytope.adjacent_facets(vertices, F, dir=test_dir)

        @test length(adj_facets) == 48

        adj_games = map(f -> convert(BellGame, f, scenario, rep="normalized"), adj_facets)
        canonical_games = unique(map(g -> LocalPolytope.generator_facet(g, scenario), adj_games))

        @test canonical_games == [
            [1 0 0 0;1 0 0 0;1 0 0 0;0 0 0 0],
            [1 0 0 0;1 0 0 0;0 1 0 0;0 0 1 0],
            [2 0 0 0;1 1 1 0;0 2 0 0;0 0 1 1],
            [1 1 0 0;1 0 1 0;0 1 1 0;0 0 0 1],
            [2 0 0 0;1 1 1 0;0 2 0 0;0 0 2 0]
        ]
    end
end

@testset "adjacency_decomposition" begin
    @testset  "41-2-14 polytope skip positivity" begin
        scenario = LocalSignaling(4,4,2)
        vertices = LocalPolytope.vertices(scenario)

        BG = BellGame([1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1], 2)

        skip = [BellGame([1 0 0 0;1 0 0 0;1 0 0 0;0 0 0 0],1)]
        dict = LocalPolytope.adjacency_decomposition(vertices, BG, scenario, skip_games=skip, dir=test_dir)

        # positvity included, but is skipped in computation
        @test collect(keys(dict)) == [
            [1 0 0 0;1 0 0 0;0 1 0 0;0 0 1 0],
            [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1],
            [1 0 0 0;1 0 0 0;1 0 0 0;0 0 0 0],
            [2 0 0 0;1 1 1 0;0 2 0 0;0 0 1 1],
            [2 0 0 0;1 1 1 0;0 2 0 0;0 0 2 0],
            [1 1 0 0;1 0 1 0;0 1 1 0;0 0 0 1],
        ]

        @test dict[[1 0 0 0;1 0 0 0;1 0 0 0;0 0 0 0]]["skipped"]

        @test all(d -> haskey(d[2], "norm_facet"), dict)
        @test all(d -> haskey(d[2], "generator_facet"), dict)
    end

    @testset  "41-2-14 polytope no skips" begin
        scenario = LocalSignaling(4,4,2)
        vertices = LocalPolytope.vertices(scenario)

        BG = BellGame([1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1], 2)

        dict = LocalPolytope.adjacency_decomposition(vertices, BG, scenario, dir=test_dir)

        @test collect(keys(dict)) == [
            [1 0 0 0;1 0 0 0;0 1 0 0;0 0 1 0],
            [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1],
            [1 0 0 0;1 0 0 0;1 0 0 0;0 0 0 0],
            [2 0 0 0;1 1 1 0;0 2 0 0;0 0 1 1],
            [2 0 0 0;1 1 1 0;0 2 0 0;0 0 2 0],
            [1 1 0 0;1 0 1 0;0 1 1 0;0 0 0 1]
        ]

        @test dict[[1 0 0 0;1 0 0 0;1 0 0 0;0 0 0 0]]["skipped"] == false
    end

    @testset "41-2-14 polytope filter out positivity vertices" begin
        scenario = LocalSignaling(4,4,2)
        vertices = LocalPolytope.vertices(scenario, rank_d_only = true)

        BG = BellGame([1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1], 2)

        dict = LocalPolytope.adjacency_decomposition(vertices, BG, scenario, dir=test_dir)

        @test issetequal( collect(keys(dict)), [
            [1 0 0 0;1 0 0 0;0 1 0 0;0 0 1 0],
            [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1],
            [1 0 0 0;1 0 0 0;1 0 0 0;0 0 0 0],
            [2 0 0 0;1 1 1 0;0 2 0 0;0 0 1 1],
            [2 0 0 0;1 1 1 0;0 2 0 0;0 0 2 0],
            [1 1 0 0;1 0 1 0;0 1 1 0;0 0 0 1],
            [1 1 1 1;0 0 0 0;0 0 0 0;0 0 0 0],
        ])
    end

    @testset "31-2-13 polytope skip by  max_vertices" begin
        scenario = LocalSignaling(3,3,2)
        vertices = LocalPolytope.vertices(scenario)

        BG  = BellGame([1 0 0;1 0 0;0 0 0], 1)

        dict = LocalPolytope.adjacency_decomposition(vertices, BG, scenario, max_vertices=5, dir=test_dir)

        @test dict[[1 0 0;0 1 0;0 0 1]]["skipped"]

        @test collect(keys(dict)) == [[1 0 0;1 0 0;0 0 0],[1 0 0;0 1 0;0 0 1]]
    end

    @testset "json logging" begin
        filename = "JSON_logging_test.json"
        if isfile(test_dir*filename)
            rm(test_dir*filename, force=true)
        end
        println("past file check")

        scenario = LocalSignaling(3,3,2)
        vertices = LocalPolytope.vertices(scenario)

        BG  = BellGame([1 0 0;1 0 0;0 0 0], 1)

        println("logging")
        dict = LocalPolytope.adjacency_decomposition(
            vertices, BG, scenario,
            dir=test_dir, log=true, log_filename = filename
        )

        println("about to parse")
        log_json = JSON.parsefile(test_dir*filename)

        @test log_json == Dict{String,Any}(
            "[1 0 0; 0 1 0; 0 0 1]" => Dict{String,Any}(
                "generator_facet" => Any[1, 0, 0, 1, -1, -1, 1],
                "norm_facet" => Any[1, 0, -1, -1, 0, 1, 1],
                "considered" => true,
                "num_vertices" => 6,
                "skipped" => false
            ),
            "[1 0 0; 1 0 0; 0 0 0]" => Dict{String,Any}(
                "generator_facet" => Any[1, 1, 0, 0, 0, 0, 1],
                "norm_facet" => Any[1, 1, 0, 0, 0, 0, 1],
                "considered" => true,
                "num_vertices" => 14,
                "skipped" => false
            )
        )

        println("about to delet file")
        rm(test_dir*filename, force=true)
    end
end

end
