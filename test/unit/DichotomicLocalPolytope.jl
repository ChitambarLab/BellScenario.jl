using Test

@testset "/src/BellScenario/DichotomicLocalPolytope.jl" begin

using BellScenario: DichotomicLocalPolytope

@testset "bipartite_no_signaling_strategies()" begin

    @testset "16 strategies are returned" begin
        s = DichotomicLocalPolytope.bipartite_no_signaling_strategies()
        @test size(s) == (16,)

        s = DichotomicLocalPolytope.bipartite_no_signaling_strategies(subspace="p12")
        @test size(s) == (16,)
    end

    @testset "strategies are 4x4" begin
        s = DichotomicLocalPolytope.bipartite_no_signaling_strategies()
        for i in 1:size(s)[1]
            @test size(s[i]) == (4,4)
        end
    end

    @testset "strategies are 3x4 in p12 subspace" begin
        s = DichotomicLocalPolytope.bipartite_no_signaling_strategies(subspace="p12")
        for i in 1:size(s)[1]
            @test size(s[i]) == (3,4)
        end
    end

    @testset "returned strategies are unique" begin
        s = DichotomicLocalPolytope.bipartite_no_signaling_strategies()
        for i in 1:size(s)[1]
            strategy = pop!(s)
            non_unique_count = 0
            for j in 1:size(s)[1]
                if strategy == s[j]
                    non_unique_count += 1
                end
            end
            @test non_unique_count == 0
        end
    end

    @testset "returned p12 strategies are unique" begin
        s = DichotomicLocalPolytope.bipartite_no_signaling_strategies(subspace="p12")
        for i in 1:size(s)[1]
            strategy = pop!(s)
            non_unique_count = 0
            for j in 1:size(s)[1]
                if strategy == s[j]
                    non_unique_count += 1
                end
            end
            @test non_unique_count == 0
        end
    end

    @testset "spot checking returned values" begin
        s = DichotomicLocalPolytope.bipartite_no_signaling_strategies()
        @test s[1] == [
            1 0 0 0;
            0 1 0 0;
            0 0 1 0;
            0 0 0 1;
        ]
        @test s[2] == [
            1 1 0 0;
            0 0 0 0;
            0 0 1 1;
            0 0 0 0;
        ]
        @test s[15] == [
            0 0 0 0;
            0 0 1 1;
            0 0 0 0;
            1 1 0 0;
        ]
    end

    @testset "spot checking returned values for p12" begin
        s = DichotomicLocalPolytope.bipartite_no_signaling_strategies(subspace="p12")
        @test s[1] == [
            1 0 0 0;
            0 1 0 0;
            0 0 1 0;
        ]
        @test s[2] == [
            1 1 0 0;
            0 0 0 0;
            0 0 1 1;
        ]
        @test s[15] == [
            0 0 0 0;
            0 0 1 1;
            0 0 0 0;
        ]
    end
end


@testset "deterministic_fixed_direction_signaling_behaviors()" begin
    test_cases = [
        ("p16",16), # (subspace, dimension of returned behaviors)
        ("p12",12),
        ("p10",10)
    ]

    @testset "subspace=$(case[1])" for case in test_cases
        vertices = DichotomicLocalPolytope.deterministic_fixed_direction_signaling_behaviors(subspace=case[1])

        @test size(vertices) == (64,)

        @testset "behaviors are unique" for v_test in vertices
            @test size(v_test) == (case[2],)

            duplicate_count = 0
            for v_match in vertices
                if v_test == v_match
                    duplicate_count += 1
                end
            end
            @test duplicate_count == 1
        end
    end

    @testset "verify equivalence between subspaces" begin
        p16_vertices = DichotomicLocalPolytope.deterministic_fixed_direction_signaling_behaviors(subspace="p16")
        p12_vertices = DichotomicLocalPolytope.deterministic_fixed_direction_signaling_behaviors(subspace="p12")
        p10_vertices = DichotomicLocalPolytope.deterministic_fixed_direction_signaling_behaviors(subspace="p10")
        for i in 1:64
            catv_16 = cat([1], p16_vertices[i], dims = 1)
            catv_12 = cat([1], p12_vertices[i], dims = 1)
            catv_10 = cat([1], p10_vertices[i], dims = 1)

            @test catv_10 == Behavior.gen_to_fd_proj((2,2),(2,2))*catv_16
            @test catv_12 == Behavior.gen_to_norm_proj((2,2),(2,2))*catv_16
            @test catv_16 == Behavior.fd_to_gen_proj((2,2),(2,2))*catv_10
            @test catv_16 == Behavior.norm_to_gen_proj((2,2),(2,2))*catv_12
        end
    end

    @testset "verify that p16_vertices are normalized" begin
        p16_vertices = DichotomicLocalPolytope.deterministic_fixed_direction_signaling_behaviors(subspace="p16")

        for v in p16_vertices
            for i in 1:4
                @test 1 == v[i] + v[i+4] + v[i+8] + v[i+12]
            end
        end
    end
end

@testset "deterministic_one_way_signaling_behaviors()" begin
    vertices = DichotomicLocalPolytope.deterministic_one_way_signaling_behaviors()

    @test size(vertices) == (112,)

    p12_vertices = DichotomicLocalPolytope.deterministic_one_way_signaling_behaviors(subspace="p12")
    @test size(p12_vertices) == (112,)

    @testset "behaviors are normalized" begin
        for vertex in vertices
            @test size(vertex) == (16,)
            # there are 4 normalization constraints per behavior
            for i in 1:4
                @test 1 == vertex[i] + vertex[4+i] + vertex[8+i] + vertex[12+i]
            end
        end

        for vertex in p12_vertices
            @test size(vertex) == (12,)
            for i in 1:4
                @test 1 >= vertex[i] + vertex[4+i] + vertex[8+i] >= 0
            end
        end
    end

    @testset "behaviors are unique" begin
        for test_vertex in vertices
            duplicate_count = 0
            for match_vertex in vertices
                if test_vertex == match_vertex
                    duplicate_count += 1
                end
            end
            @test duplicate_count == 1
        end

        for test_vertex in p12_vertices
            duplicate_count = 0
            for match_vertex in p12_vertices
                if test_vertex == match_vertex
                    duplicate_count += 1
                end
            end
            @test duplicate_count == 1
        end
    end
end

@testset "deterministic_no_signaling_behaviors()" begin

    @testset "16 vertices are returned for bipartite case" begin
        vertices = DichotomicLocalPolytope.deterministic_no_signaling_behaviors()
        @test size(vertices) == (16,)

        p12_vertices = DichotomicLocalPolytope.deterministic_no_signaling_behaviors(subspace="p12")
        @test size(p12_vertices) == (16,)
    end

    @testset "vertices are normalized with respect to unique inputs" begin
        vertices = DichotomicLocalPolytope.deterministic_no_signaling_behaviors()
        for vertex in vertices
            # there are 4 normalization constraints per behavior
            for i in 1:4
                @test 1 == vertex[i] + vertex[4+i] + vertex[8+i] + vertex[12+i]
            end
        end

        p12_vertices = DichotomicLocalPolytope.deterministic_no_signaling_behaviors(subspace="p12")
        for vertex in p12_vertices
            for i in 1:4
                @test 1 >= vertex[i] + vertex[4+i] + vertex[8+i] >= 0
            end
        end
    end

    @testset "vertices are unique" begin
        vertices = DichotomicLocalPolytope.deterministic_no_signaling_behaviors()

        for i in 1:size(vertices)[1]
            test_vertex = pop!(vertices)
            non_unique_count = 0
            for vertex in vertices
                if test_vertex == vertex
                    non_unique_count += 1
                end
            end
            @test non_unique_count == 0
        end

        p12_vertices = DichotomicLocalPolytope.deterministic_no_signaling_behaviors(subspace="p12")

        for test_vertex in p12_vertices
            duplicate_count = 0
            for match_vertex in p12_vertices
                if test_vertex == match_vertex
                    duplicate_count += 1
                end
            end
            @test duplicate_count == 1
        end
    end

    @testset "p16 vertice match p12 vertices" begin
        p16_vertices = DichotomicLocalPolytope.deterministic_no_signaling_behaviors(subspace="p16")
        p12_vertices = DichotomicLocalPolytope.deterministic_no_signaling_behaviors(subspace="p12")

        p12_match = []
        for vertex in p16_vertices
            catv_16 = cat([1],vertex,dims=1)
            push!(p12_match, (Behavior.gen_to_norm_proj((2,2),(2,2))*catv_16)[2:end] )
        end
        @test p12_vertices == p12_match

        p16_match = []
        for vertex in p12_vertices
            catv_12 = cat([1],vertex,dims=1)
            push!(p16_match, (Behavior.norm_to_gen_proj((2,2),(2,2))*catv_12)[2:end])
        end
        @test p16_vertices == p16_match
    end

    @testset "vertices match the no-signaling polytope" begin
        vertices = DichotomicLocalPolytope.deterministic_no_signaling_behaviors()

        no_signaling_vertices_match = [
            [0,0,0,0,0,0,0,0],
            [0,0,0,1,0,0,0,0],
            [0,0,1,0,0,0,0,0],
            [0,0,1,1,0,0,0,0],
            [0,1,0,0,0,0,0,0],
            [0,1,0,1,0,0,0,1],
            [0,1,1,0,0,0,1,0],
            [0,1,1,1,0,0,1,1],
            [1,0,0,0,0,0,0,0],
            [1,0,0,1,0,1,0,0],
            [1,0,1,0,1,0,0,0],
            [1,0,1,1,1,1,0,0],
            [1,1,0,0,0,0,0,0],
            [1,1,0,1,0,1,0,1],
            [1,1,1,0,1,0,1,0],
            [1,1,1,1,1,1,1,1],
        ]

        # loop through each computed vertex and make sure that it is in the
        # match set only once.
        for vertex in vertices
            catv = cat([1],vertex,dims=1)
            no_signaling_vertex = (Behavior.gen_to_ns_proj((2,2),(2,2))*catv)[2:end]

            match_found = 0
            for match_vertex in no_signaling_vertices_match
                if no_signaling_vertex == match_vertex
                    match_found += 1
                end
            end
            @test match_found == 1
        end
    end
end

@testset "two_input_strategies()" begin
    strategies = DichotomicLocalPolytope.two_input_strategies()

    @test size(strategies) == (16,)

    @testset "two_input_strategies are unique amongst themselves" begin
        for test_strategy in strategies
            duplicate_count = 0
            for match_strategy in strategies
                if test_strategy == match_strategy
                    duplicate_count += 1
                end
            end
            # test_strategy should only match with itself
            @test duplicate_count == 1
        end
    end
end

@testset "dichotomic_strategies()" begin
    s = DichotomicLocalPolytope.dichotomic_strategies()
    @test s[1] == [1 0; 0 1]
    @test s[2] == [1 1; 0 0]
    @test s[3] == [0 0; 1 1]
    @test s[4] == [0 1; 1 0]
end

@testset "Permutations of communication protocols create redundant strategies" begin
    # consider a two-input strategy α(x,σ(y)) → α*(x⊗σ(y)) = α*(I⊗σ)*(x⊗y)
    # in the below tests verify that { α } = { α*(I⊗σ) } i.e. the set of
    # permuted strategies contains all members of {α} and that {α} contains the
    # set of all permuted strategies.

    two_input_strategies = DichotomicLocalPolytope.two_input_strategies()
    dichotomic_strategies = DichotomicLocalPolytope.dichotomic_strategies()

    permuted_strategies = []
    for ds in dichotomic_strategies[[2 3 4]]
        input_matrix = kron([1 0; 0 1], ds)
        for ts in two_input_strategies
            push!(permuted_strategies,ts*input_matrix)
        end
    end

    @testset "each permuted strategy appears once in the set of two input strategies" begin
        for ps_test in permuted_strategies
            duplicate_count = 0
            for ts_match in two_input_strategies
                if ps_test == ts_match
                    duplicate_count += 1
                end
            end

            @test duplicate_count == 1
        end
    end

    @testset "each two input strategy appears at least once in the set of permuted strategies" begin
        for ts_test in two_input_strategies
            duplicate_count = 0
            for ps_match in permuted_strategies
                if ts_test == ps_match
                    duplicate_count += 1
                end
            end
            @test duplicate_count > 0
        end
    end
end

end
