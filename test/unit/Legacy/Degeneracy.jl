using Test, LinearAlgebra
using QBase: permutation_matrices

@testset "/src/Legacy/Degeneracy.jl" begin

using BellScenario: Degeneracy

@testset "symmetry_groups()" begin
    @testset "2x2 permutation group" begin
        transforms = [[1 0;0 1],[0 1;1 0]]
        elements = [[0 0],[1 0],[0 1],[1 1]]

        sym_groups = Degeneracy.symmetry_groups(transforms,elements)

        @test issetequal(elements, collect(Iterators.flatten(sym_groups)))
        @test length(sym_groups) == 3

        @test sym_groups[1] == [[0 0]]
        @test sym_groups[2] == [[1 0],[0 1]]
        @test sym_groups[3] == [[1 1]]
    end

    @testset "3x3 permutation group" begin
        transforms = permutation_matrices(3)
        elements = [
            [0 0 0],[0 0 1],[0 1 0],[0 1 1],
            [1 0 0],[1 0 1],[1 1 0],[1 1 1]
        ]

        sym_groups = Degeneracy.symmetry_groups(transforms,elements)

        @test issetequal(elements, collect(Iterators.flatten(sym_groups)))
        @test length(sym_groups) == 4

        @test sym_groups[1] == [[0 0 0]]
        @test sym_groups[2] == [[0 0 1],[0 1 0],[1 0 0]]
        @test sym_groups[3] == [[0 1 1],[1 0 1],[1 1 0]]
        @test sym_groups[4] == [[1 1 1]]
    end

    @testset "4x4 permutation group" begin
        transforms = permutation_matrices(4)
        elements = [
            [0 0 0 0],[0 0 0 1],[0 0 1 0],[0 0 1 1],
            [0 1 0 0],[0 1 0 1],[0 1 1 0],[0 1 1 1],
            [1 0 0 0],[1 0 0 1],[1 0 1 0],[1 0 1 1],
            [1 1 0 0],[1 1 0 1],[1 1 1 0],[1 1 1 1]
        ]

        sym_groups = Degeneracy.symmetry_groups(transforms,elements)

        @test issetequal(elements, collect(Iterators.flatten(sym_groups)))
        @test length(sym_groups) == 5

        @test sym_groups[1] == [[0 0 0 0]]
        @test sym_groups[2] == [[0 0 0 1],[0 0 1 0],[0 1 0 0],[1 0 0 0]]
        @test sym_groups[3] == [[0 0 1 1],[0 1 0 1],[0 1 1 0],[1 0 0 1],[1 0 1 0],[1 1 0 0]]
        @test sym_groups[4] == [[0 1 1 1],[1 0 1 1],[1 1 0 1],[1 1 1 0]]
        @test sym_groups[5] == [[1 1 1 1]]
    end
end

@testset "merge_symmetry_groups()" begin
    group1_sets = [[[1,1],[2,2]],[[3,3],[4,4],[5,5]],[[6,6]],[[7,7]]]
    group2_sets = [[[1,1],[3,3]],[[6,6],[7,7]],[[2,2],[4,4],[5,5]]]

    merged_groups = Degeneracy.merge_symmetry_groups(group1_sets,group2_sets)
    @test length(merged_groups) == 2
    @test all( # returned group sets do not overlap
        (x) -> x === nothing,
        indexin(merged_groups[1],merged_groups[2])
    )
    @test issetequal(merged_groups[1], [[1,1],[2,2],[3,3],[4,4],[5,5]])
    @test issetequal(merged_groups[2], [[6,6],[7,7]])
end

@testset "canonical_facets()" begin
    α_expt = (2,2)
    β_expt = (2,2)

    facets = [
        [0 0 0 0 0 -1 0 0 0],[0 0 0 0 0 0 -1 0 0],[0 0 0 0 0 0 0 -1 0],
        [0 0 0 0 0 0 0 0 -1],[0 -1 0 0 0 0 1 0 0],[0 -1 0 0 0 1 0 0 0],
        [0 0 -1 0 0 0 0 0 1],[0 0 -1 0 0 0 0 1 0],[0 0 0 -1 0 0 0 1 0],
        [0 0 0 -1 0 1 0 0 0],[0 0 0 0 -1 0 0 0 1],[0 0 0 0 -1 0 1 0 0],
        [0 -1 0 -1 0 1 1 1 -1],[0 -1 0 0 -1 1 1 -1 1],[0 0 -1 -1 0 1 -1 1 1],
        [0 0 -1 0 -1 -1 1 1 1],[-1 0 1 0 1 0 0 0 -1],[-1 0 1 1 0 0 0 -1 0],
        [-1 1 0 0 1 0 -1 0 0],[-1 1 0 1 0 -1 0 0 0],[-1 0 1 0 1 1 -1 -1 -1],
        [-1 0 1 1 0 -1 1 -1 -1],[-1 1 0 0 1 -1 -1 1 -1],[-1 1 0 1 0 -1 -1 -1 1]
    ]

    canonical_groups = Degeneracy.canonical_facets(α_expt, β_expt, facets, "non-signaling")

    @test issetequal(facets, collect(Iterators.flatten(canonical_groups)))
    @test length(canonical_groups) == 2
    @test issetequal(canonical_groups[1], [
        [0 0 0 0 0 -1 0 0 0],[0 0 0 0 0 0 -1 0 0],[0 0 0 0 0 0 0 -1 0],[0 0 0 0 0 0 0 0 -1],
        [0 -1 0 0 0 0 1 0 0],[0 -1 0 0 0 1 0 0 0],[0 0 -1 0 0 0 0 0 1],[0 0 -1 0 0 0 0 1 0],
        [0 0 0 -1 0 0 0 1 0],[0 0 0 -1 0 1 0 0 0],[0 0 0 0 -1 0 0 0 1],[0 0 0 0 -1 0 1 0 0],
        [-1 0 1 0 1 0 0 0 -1],[-1 0 1 1 0 0 0 -1 0],[-1 1 0 0 1 0 -1 0 0],[-1 1 0 1 0 -1 0 0 0]
    ])
    @test issetequal(canonical_groups[2], [
        [0 -1 0 -1 0 1 1 1 -1],[0 -1 0 0 -1 1 1 -1 1],[0 0 -1 -1 0 1 -1 1 1],[0 0 -1 0 -1 -1 1 1 1],
        [-1 0 1 0 1 1 -1 -1 -1],[-1 0 1 1 0 -1 1 -1 -1],[-1 1 0 0 1 -1 -1 1 -1],[-1 1 0 1 0 -1 -1 -1 1]
    ])
end

@testset "canonical_vertices()" begin
    @testset "chsh vertices" begin
        α1 = (2,2)
        β1 = (2,2)

        vertices = [
            [1;0;0;0;0;0;0;0;0],[1;0;0;0;1;0;0;0;0],[1;0;0;1;0;0;0;0;0],[1;0;0;1;1;0;0;0;0],
            [1;0;1;0;0;0;0;0;0],[1;0;1;0;1;0;0;0;1],[1;0;1;1;0;0;0;1;0],[1;0;1;1;1;0;0;1;1],
            [1;1;0;0;0;0;0;0;0],[1;1;0;0;1;0;1;0;0],[1;1;0;1;0;1;0;0;0],[1;1;0;1;1;1;1;0;0],
            [1;1;1;0;0;0;0;0;0],[1;1;1;0;1;0;1;0;1],[1;1;1;1;0;1;0;1;0],[1;1;1;1;1;1;1;1;1]
        ]

        canonical_vertices = Degeneracy.canonical_vertices(α1, β1, vertices, "non-signaling")

        @test length(canonical_vertices) == 1
        @test length(canonical_vertices[1]) == 16
        @test issetequal(canonical_vertices[1], vertices)
    end

    @testset "prepare and measure vertices" begin
        α = (3,1)
        β = (1,3)

        vertices = [
            [1 1 0 0 0 0 0]',[1 1 1 0 0 0 0]',[1 1 0 1 0 0 0]',[1 1 1 1 0 0 0]',[1 1 1 0 0 0 1]',[1 1 0 1 0 1 0]',
            [1 1 0 0 0 1 1]',[1 0 0 0 1 0 0]',[1 0 0 0 1 1 0]',[1 0 0 0 1 0 1]',[1 0 1 1 1 0 0]',[1 0 1 0 1 0 1]',
            [1 0 0 1 1 1 0]',[1 0 0 0 1 1 1]',[1 0 0 0 0 0 0]',[1 0 1 0 0 0 0]',[1 0 0 1 0 0 0]',[1 0 1 1 0 0 0]',
            [1 0 0 0 0 1 1]',[1 0 0 0 0 1 0]',[1 0 0 0 0 0 1]'
        ]

        canonical_vertices = Degeneracy.canonical_vertices(α, β, vertices, "normalized")

        @test length(canonical_vertices) == 2
        @test length(canonical_vertices[1]) == 3
        @test length(canonical_vertices[2]) == 18
        @test issetequal(canonical_vertices[1], [[1 1 1 1 0 0 0]',[1 0 0 0 0 0 0]',[1 0 0 0 1 1 1]'])
        @test issetequal( canonical_vertices[2],[
            [1 1 0 0 0 0 0]',[1 1 1 0 0 0 0]',[1 1 0 1 0 0 0]',[1 1 1 0 0 0 1]',[1 1 0 1 0 1 0]',
            [1 1 0 0 0 1 1]',[1 0 0 0 1 0 0]',[1 0 0 0 1 1 0]',[1 0 0 0 1 0 1]',[1 0 1 1 1 0 0]',
            [1 0 0 1 1 1 0]',[1 0 1 0 0 0 0]',[1 0 0 1 0 0 0]',[1 0 1 1 0 0 0]',[1 0 1 0 1 0 1]',
            [1 0 0 0 0 1 1]',[1 0 0 0 0 1 0]',[1 0 0 0 0 0 1]'
        ])
    end
end

@testset "bipartite_input_relabels()" begin

    @testset "non-signaling subspace" begin
        id = diagm(0 => fill(1, 9))
        flip_alice_input = [
            1 0 0 0 0 0 0 0 0;
            0 0 1 0 0 0 0 0 0;
            0 1 0 0 0 0 0 0 0;
            0 0 0 1 0 0 0 0 0;
            0 0 0 0 1 0 0 0 0;
            0 0 0 0 0 0 0 1 0;
            0 0 0 0 0 0 0 0 1;
            0 0 0 0 0 1 0 0 0;
            0 0 0 0 0 0 1 0 0;
        ];
        flip_bob_input = [
            1 0 0 0 0 0 0 0 0;
            0 1 0 0 0 0 0 0 0;
            0 0 1 0 0 0 0 0 0;
            0 0 0 0 1 0 0 0 0;
            0 0 0 1 0 0 0 0 0;
            0 0 0 0 0 0 1 0 0;
            0 0 0 0 0 1 0 0 0;
            0 0 0 0 0 0 0 0 1;
            0 0 0 0 0 0 0 1 0;
        ];
        flip_ab_input = [
            1 0 0 0 0 0 0 0 0;
            0 0 1 0 0 0 0 0 0;
            0 1 0 0 0 0 0 0 0;
            0 0 0 0 1 0 0 0 0;
            0 0 0 1 0 0 0 0 0;
            0 0 0 0 0 0 0 0 1;
            0 0 0 0 0 0 0 1 0;
            0 0 0 0 0 0 1 0 0;
            0 0 0 0 0 1 0 0 0;
        ]

        @test issetequal(
            Degeneracy.bipartite_input_relabels((2,2), (2,2), "non-signaling"),
            [id, flip_alice_input, flip_bob_input, flip_ab_input]
        )
    end

    test_cases = [
        ("non-signaling", 8),
        ("fixed-direction", 10),
        ("normalized", 12),
        ("generalized", 16)
    ]
    @testset "basic tests for subspace: $(case[1])" for case in test_cases
        maps = Degeneracy.bipartite_input_relabels((2,2),(2,2),case[1])

        dim = case[2]

        @test size(maps) == size(unique(maps))

        for i in 1:size(maps)[1]
            map = maps[i]

            @test size(map) == (dim+1,dim+1)
            @test map[1,1] == 1
            @test map[2:end,1] == zeros(Int64, dim)
            @test map[1,2:end] == zeros(Int64, dim)
        end
    end
end

@testset "bipartite_output_relabels()" begin
    @testset "dichotomic case" begin
        @testset "non-signaling" begin
            (α_relabels, β_relabels) = Degeneracy.bipartite_output_relabels((2,2),(2,2),"non-signaling")

            @test size(α_relabels) == (4,)
            @test size(β_relabels) == (4,)

            @test α_relabels[1] == diagm(0=>fill(1, 9))
            @test β_relabels[1] == diagm(0=>fill(1, 9))

            @test α_relabels[2] == [
                1 0 0 0 0 0  0 0 0;
                1 -1 0 0 0 0 0 0 0;
                0 0 1 0 0 0 0 0 0;
                0 0 0 1 0 0 0 0 0;
                0 0 0 0 1 0 0 0 0;
                0 0 0 1 0 -1 0 0 0;
                0 0 0 0 1 0 -1 0 0;
                0 0 0 0 0 0 0 1 0;
                0 0 0 0 0 0 0 0 1;
            ]
            @test α_relabels[3] == [
                 1 0 0 0 0 0 0 0 0;
                 0 1 0 0 0 0 0 0 0;
                 1 0 -1 0 0 0 0 0 0;
                 0 0 0 1 0 0 0 0 0;
                 0 0 0 0 1 0 0 0 0;
                 0 0 0 0 0 1 0 0 0;
                 0 0 0 0 0 0 1 0 0;
                 0 0 0 1 0 0 0 -1 0;
                 0 0 0 0 1 0 0 0 -1;
             ]
             @test β_relabels[2] == [
                1 0 0 0 0 0 0 0 0;
                0 1 0 0 0 0 0 0 0;
                0 0 1 0 0 0 0 0 0;
                1 0 0 -1 0 0 0 0 0;
                0 0 0 0 1 0 0 0 0;
                0 1 0 0 0 -1 0 0 0;
                0 0 0 0 0 0 1 0 0;
                0 0 1 0 0 0 0 -1 0;
                0 0 0 0 0 0 0 0 1;
            ]
            @test β_relabels[3] == [
                1 0 0 0 0 0 0 0 0;
                0 1 0 0 0 0 0 0 0;
                0 0 1 0 0 0 0 0 0;
                0 0 0 1 0 0 0 0 0;
                1 0 0 0 -1 0 0 0 0;
                0 0 0 0 0 1 0 0 0;
                0 1 0 0 0 0 -1 0 0;
                0 0 0 0 0 0 0 1 0;
                0 0 1 0 0 0 0 0 -1;
            ]
            @test β_relabels[4] == β_relabels[2]*β_relabels[3]
        end
    end

    @testset "relabel set contains the inverse of each element" begin
        test_cases = [
            ((2,2),(2,2)),
            ((3,2),(2,2)),
            ((2,3),(2,3)),
            ((3,1),(1,3))
        ]

        @testset "generalized $(case[1]), $(case[2])" for case in test_cases
            (α_relabels, β_relabels) = Degeneracy.bipartite_output_relabels(case[1],case[2],"generalized")

            dim = Behavior.dimension(case[1],case[2],"generalized")

            # in generalized coordinates, output relabelings inverted through
            # their own transpose.
            for relabel in α_relabels
                parity = det(relabel)
                @test (parity == 1) | (parity == -1)
                @test relabel*relabel' == diagm( 0 => fill(1,dim))
            end

            for relabel in β_relabels
                parity = det(relabel)
                @test (parity == 1) | (parity == -1)
                @test relabel*relabel' == diagm( 0 => fill(1,dim))
            end
        end

        @testset "non-signaling $(case[1]), $(case[2])" for case in test_cases
            (α_relabels, β_relabels) = Degeneracy.bipartite_output_relabels(case[1],case[2],"non-signaling")

            dim = Behavior.dimension(case[1],case[2],"non-signaling")
            Id = diagm(0 => fill(1,dim))

            foreach( relabel_test->begin
                matches = filter(relabel_inv->relabel_test*relabel_inv == Id, α_relabels)
                @test (1,) == size(matches)
            end, α_relabels)

            foreach( relabel_test->begin
                matches = filter(relabel_inv->relabel_test*relabel_inv == Id, β_relabels)
                @test (1,) == size(matches)
            end, β_relabels)
        end

        @testset "fixed-direction $(case[1]), $(case[2])" for case in test_cases
            (α_relabels, β_relabels) = Degeneracy.bipartite_output_relabels(case[1],case[2],"fixed-direction")

            dim = Behavior.dimension(case[1],case[2],"fixed-direction")
            Id = diagm(0 => fill(1,dim))

            # check that there is exactly 1 inverting map for every map in the
            # relabeling set
            foreach( relabel_test->begin
                matches = filter(relabel_inv->relabel_test*relabel_inv == Id, α_relabels)
                @test (1,) == size(matches)
            end, α_relabels)

            foreach( relabel_test->begin
                matches = filter(relabel_inv->relabel_test*relabel_inv == Id, β_relabels)
                @test (1,) == size(matches)
            end, β_relabels)
        end
    end
end

end
