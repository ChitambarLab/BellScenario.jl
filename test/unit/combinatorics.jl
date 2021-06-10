using Test, LinearAlgebra
using BellScenario: stirling2, stirling2_partitions, stirling2_matrices, permutation_matrices, base_n_val, n_choose_k_matrices

@testset "./src/combinatorics.jl" begin

@testset "stirling2()" begin
    @testset "domain errors" begin
        @test_throws DomainError stirling2(1,0)
        @test_throws DomainError stirling2(1,2)
    end

    @test stirling2(1,1) == 1

    @testset "n = 2" begin
        @test stirling2(2,1) == 1
        @test stirling2(2,2) == 1
    end

    @testset "n = 3" begin
        @test stirling2(3,1) == 1
        @test stirling2(3,2) == 3
        @test stirling2(3,3) == 1
    end

    @testset "n = 4" begin
        @test stirling2(4,1) == 1
        @test stirling2(4,2) == 7
        @test stirling2(4,3) == 6
        @test stirling2(4,4) == 1
    end

    @testset "n = 5" begin
        @test stirling2(5,1) == 1
        @test stirling2(5,2) == 15
        @test stirling2(5,3) == 25
        @test stirling2(5,4) == 10
        @test stirling2(5,5) == 1
    end

    @testset "n = 6" begin
        @test stirling2(6,1) == 1
        @test stirling2(6,2) == 31
        @test stirling2(6,3) == 90
        @test stirling2(6,4) == 65
        @test stirling2(6,5) == 15
        @test stirling2(6,6) == 1
    end

    @testset "n = 10" begin
        @test stirling2(10,1) == 1
        @test stirling2(10,2) == 511
        @test stirling2(10,3) == 9330
        @test stirling2(10,4) == 34105
        @test stirling2(10,5) == 42525
        @test stirling2(10,6) == 22827
        @test stirling2(10,7) == 5880
        @test stirling2(10,8) == 750
        @test stirling2(10,9) == 45
        @test stirling2(10,10) == 1
    end
end

@testset "stirling2_partitions()" begin
    @testset "base cases" begin
        @test stirling2_partitions(0,0) == []
        @test stirling2_partitions(1,0) == []
        @test stirling2_partitions(1,2) == []
        @test stirling2_partitions(5,1) == [[[1,2,3,4,5]]]
    end

    @testset "simple cases" begin
        @test stirling2_partitions(3,2) == [[[1,2],[3]],[[2],[1,3]],[[1],[2,3]]]

        @test stirling2_partitions(4,2) == [
            [[1, 2, 3], [4]],
            [[3], [1, 2, 4]],
            [[1, 2], [3, 4]],
            [[1, 3], [2, 4]],
            [[2], [1, 3, 4]],
            [[2, 3], [1, 4]],
            [[1], [2, 3, 4]],
        ]

        @test stirling2_partitions(4,4) == [[[1],[2],[3],[4]]]
    end

    @testset "general partitions" begin
        for n in 2:7
            for k in 1:n
                partitions = stirling2_partitions(n,k)
                @test length(partitions) == stirling2(n,k)
                @test length(unique(partitions)) == stirling2(n,k)
                @test all(s -> length(s) == k, partitions)
                @test all(s -> sort(cat(s...,dims=1)) == [1:n...], partitions)
            end
        end
    end
end

@testset "stirling2_matrices()" begin
    @testset "domain errors" begin
        @test_throws DomainError n_choose_k_matrices(2,3)
        @test_throws DomainError n_choose_k_matrices(2,0)
    end

    @testset "simple cases" begin
        @test stirling2_matrices(4,1) == [
            [1 1 1 1]
        ]
        @test stirling2_matrices(3,2) == [
            [1 1 0;0 0 1],[0 1 0;1 0 1],[1 0 0;0 1 1]
        ]
    end

    @testset "general cases" begin
        for n in 2:7
            for k in 1:n
                matrices = stirling2_matrices(n,k)

                @test length(matrices) == stirling2(n,k)
                @test length(unique(matrices)) == stirling2(n,k)
            end
        end
    end
end

@testset "permutation_matrices()" begin
    @testset "dim = 2" begin
        @test permutation_matrices(2) == [[1 0; 0 1], [0 1; 1 0]]
    end

    @testset "dim = 3" begin
        @test issetequal(permutation_matrices(3), [
            [1 0 0; 0 1 0; 0 0 1],
            [1 0 0; 0 0 1; 0 1 0],
            [0 1 0; 1 0 0; 0 0 1],
            [0 1 0; 0 0 1; 1 0 0],
            [0 0 1; 1 0 0; 0 1 0],
            [0 0 1; 0 1 0; 1 0 0],
        ])
    end

    @testset "permtutation map validity" for N in 4:5
        perms = permutation_matrices(N)

        @testset "maps are unique" begin
            @test size(perms) == (factorial(N),)
            @test size(perms) == size(unique(perms))
        end

        @testset "matrix dimension, determinant, and invertibility" for perm in perms
            # dimension
            @test size(perm) == (N,N)

            # parity
            parity = det(perm)
            @test (parity == -1) || (parity == 1)

            # invertibility
            @test perm*perm' == Matrix{Int64}(I,N,N)
        end
    end
end

@testset "n_choose_k_matrices()" begin
    @testset "Domain Errors" begin
        @test_throws DomainError n_choose_k_matrices(2,3)
        @test_throws DomainError n_choose_k_matrices(2,0)
    end

    @testset "simple testcases" begin
        @test n_choose_k_matrices(3,3) == [[1 0 0;0 1 0;0 0 1]]

        @test n_choose_k_matrices(3,2) == [
            [1 0;0 1;0 0],[1 0;0 0;0 1],[0 0;1 0;0 1]
        ]

        @test n_choose_k_matrices(4,2) == [
            [1 0;0 1;0 0;0 0],[1 0;0 0;0 1;0 0],[1 0;0 0;0 0;0 1],
            [0 0;1 0;0 1;0 0],[0 0;1 0;0 0;0 1],[0 0;0 0;1 0;0 1]
        ]
    end

    @testset "general tests" begin
        for n in 2:5
            for k in 1:n
                matrices = n_choose_k_matrices(n,k)

                @test length(matrices) == binomial(n,k)
                @test length(unique(matrices)) == binomial(n,k)
            end
        end
    end
end

@testset "base_n_val()" begin
    @testset "little endian" begin
        @test 1 == base_n_val(digits(1,base=2,pad=2), 2, big_endian = false)
        @test 57 == base_n_val(digits(57,base=3,pad=5), 3, big_endian = false)
    end

    @testset "big endian" begin
        @test 12 == base_n_val([1,1,0,0],2)
        @test 974 == base_n_val([1,2,3,4,4],5)
    end

    @testset "argument errors" begin
        @test_throws ArgumentError base_n_val([1,2,0,5],5)
        @test_throws ArgumentError base_n_val([1,2,0,5],4)
        @test_throws ArgumentError base_n_val([-1,0,1],2)
    end
end

end
