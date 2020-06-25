using Test, LinearAlgebra

@testset "test/unit/LocalPolytope.jl" begin

using BellScenario: LocalPolytope, DichotomicLocalPolytope, Behavior

@testset "strategies()" begin

    @testset "trivial strategies" begin
        @test LocalPolytope.strategies(1,1) == [fill(1,(1,1))]
        @test LocalPolytope.strategies(1,2) == [[1 0]',[0 1]']
        @test LocalPolytope.strategies(2,1) == [[1 1]]
        @test LocalPolytope.strategies(3,1) == [[1 1 1]]
    end

    @testset "dichotomic strategies" begin
        s = LocalPolytope.strategies(2,2)
        @test size(s) == (4,)

        @test s[1] == [1 1; 0 0]
        @test s[2] == [0 1; 1 0]
        @test s[3] == [1 0; 0 1]
        @test s[4] == [0 0; 1 1]
    end

    @testset "two_input_strategies are unique amongst themselves" begin
        strategies = LocalPolytope.strategies(4,2)
        @test size(strategies) == (16,)
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

    @testset "checking size of set" begin
        s = LocalPolytope.strategies(10,3)
        @test size(s) == (59049,)
        @test size(s[1]) == (3,10)

        s = LocalPolytope.strategies(3,10)
        @test size(s) == (1000,)
        @test size(s[1]) == (10,3)

        s = LocalPolytope.strategies(5,5)
        @test size(s) == (3125,)
        @test size(s[1]) == (5,5)
    end
end

@testset "input_set()" begin
    LocalPolytope.input_set(1)
    @testset "num_inputs = $N" for N in 1:10
        input_set = LocalPolytope.input_set(N)
        test_matrix = zeros(N,N)
        for i in 1:N
            test_matrix[:,i] += input_set[i]
        end
        identity_matrix = diagm(0 => fill(1,N))

        # the input set should form the identity matrix
        @test test_matrix == identity_matrix
    end
end

@testset "communication_protocols()" begin

    @testset "fewer inputs than outputs" begin
        @testset "simple cases" begin
            @test LocalPolytope.communication_protocols(1,2) == [[1 0]']
            @test LocalPolytope.communication_protocols(2,3) == [[1 0;0 1;0 0]]
            @test LocalPolytope.communication_protocols(2,4) == [[1 0;0 1;0 0;0 0]]
            @test LocalPolytope.communication_protocols(3,4) == [[1 0 0;0 1 0;0 0 1;0 0 0]]
        end

        @testset "scales like binomial coefficient: outputs = $num_out " for num_out in 2:7
            @testset "inputs = $num_in" for num_in in 1:num_out
                @test length(LocalPolytope.communication_protocols(num_in,num_out)) == 1
            end
        end
    end

    @testset "trivial edge cases" begin
        @test LocalPolytope.communication_protocols(2,2) == [[1 0;0 1]]
        @test LocalPolytope.communication_protocols(3,2) == [[1 1 0;0 0 1],[1 0 1;0 1 0],[1 0 0;0 1 1]]
    end

    @testset "there exists a single protocol for $N:$N and $N:1" for N in 1:10
        @test 1 == size(LocalPolytope.communication_protocols(N,N))[1]
        @test 1 == size(LocalPolytope.communication_protocols(N,1))[1]
    end

    stirling_base_cases = [
        (3,2,3),
        (4,2,7),
        (5,2,15),
        (6,2,31),
        (7,2,63),
        (8,2,127),
        (4,3,6),
        (5,3,25),
        (6,3,90),
        (5,4,10),
        (6,4,65),
        (7,4,350),
        (6,5,15),
        (7,6,21),
        (8,7,28),
    ]

    @testset "verify stirling number base case: n=$(case[1]), k=$(case[2])" for case in stirling_base_cases
        protocols = LocalPolytope.communication_protocols(case[1],case[2])
        # protocols match stirling's number of the second kind
        @test case[3] == size(protocols)[1]
        # protocols are unique
        @test size(protocols) == size(unique(protocols))

        # protocols map a single input to a single output
        for protocol in protocols
            row_sum = zeros(case[1])
            for row_id in 1:case[2]
                row_sum += protocol[row_id,:]
            end
            @test row_sum == ones(case[1])
        end
    end
end

@testset "is_no_signaling_duplicate()" begin
    @test LocalPolytope.is_no_signaling_duplicate([1 1 1 1;0 0 0 0], 2, 2)
    @test LocalPolytope.is_no_signaling_duplicate([1 1 0 0;0 0 1 1], 2, 2)
    @test !( LocalPolytope.is_no_signaling_duplicate([1 0 1 1;0 1 0 0], 2, 2) )
    @test LocalPolytope.is_no_signaling_duplicate([1 1 1 0 0 0;0 0 0 1 1 1], 2, 3)
end

@testset "vertices()" begin
    @testset "trivial no-signaling cases" begin
        @test LocalPolytope.vertices((1,1),(1,1)) == [fill(1,(1,1))]

        @test LocalPolytope.vertices((2,1),(1,1)) == [[1 1]']
        @test LocalPolytope.vertices((1,2),(1,1)) == [fill(1,(1,1)),fill(0,(1,1))]
        @test LocalPolytope.vertices((1,1),(2,1)) == [[1 1]']
        @test LocalPolytope.vertices((1,1),(1,2)) == [fill(1,(1,1)),fill(0,(1,1))]

        @test LocalPolytope.vertices((2,2),(1,1)) == [[1 1]',[0 1]',[1 0]',[0 0]']
        @test LocalPolytope.vertices((2,1),(2,1)) == [[1 1 1 1]']
        @test LocalPolytope.vertices((2,1),(1,2)) == [[1 1]',[0 0]']
        @test LocalPolytope.vertices((1,2),(2,1)) == [[1 1]',[0 0]']
        @test LocalPolytope.vertices((1,2),(1,2)) == [
            [1 1 1]',[1 0 0]',[0 1 0]',[0 0 0]'
        ]
        @test LocalPolytope.vertices((1,1),(2,2)) == [[1 1]',[0 1]',[1 0]',[0 0]']

        @test LocalPolytope.vertices((2,2),(2,1)) == [
            [1 1 1 1]',[0 0 1 1]',[1 1 0 0]',[0 0 0 0]'
        ]
        @test LocalPolytope.vertices((2,2),(1,2)) == [
            [1 1 1 1 1]',[1 1 0 0 0]',[0 1 1 0 1]',[0 1 0 0 0]',
            [1 0 1 1 0]',[1 0 0 0 0]',[0 0 1 0 0]',[0 0 0 0 0]'
        ]
        @test LocalPolytope.vertices((2,1),(2,2)) == [
            [1 1 1 1]',[0 1 0 1]',[1 0 1 0]',[0 0 0 0]'
        ]
        @test LocalPolytope.vertices((1,2),(2,2)) == [
            [1 1 1 1 1]',[1 0 1 0 1]',[1 1 0 1 0]',[1 0 0 0 0]',
            [0 1 1 0 0]',[0 0 1 0 0]',[0 1 0 0 0]',[0 0 0 0 0]'
        ]

        @test LocalPolytope.vertices((3,1),(1,2)) == [[1 1 1]',[0 0 0]']
        @test LocalPolytope.vertices((3,1),(1,3)) == [[1 1 1 0 0 0]',[0 0 0 1 1 1]',[0 0 0 0 0 0]']
    end

    @testset "trivial fixed-direction signaling cases" begin
        @test LocalPolytope.vertices((1,1),(1,1),bits=1) == [fill(1,(1,1))]

        @test LocalPolytope.vertices((2,1),(1,1),bits=1) == [[1 1]']
        @test LocalPolytope.vertices((1,2),(1,1),bits=1) == [fill(1,(1,1)),fill(0,(1,1))]
        @test LocalPolytope.vertices((1,1),(2,1),bits=1) == [[1 1]']
        @test LocalPolytope.vertices((1,1),(1,2),bits=1) == [fill(1,(1,1)),fill(0,(1,1))]

        @test LocalPolytope.vertices((2,2),(1,1),bits=1) == [[1 1]',[0 1]',[1 0]',[0 0]']
        @test LocalPolytope.vertices((2,1),(2,1),bits=1) == [[1 1 1 1]']
        @test LocalPolytope.vertices((2,1),(1,2),bits=1) == [[1 1]',[0 0]',[0 1]',[1 0]']
        @test LocalPolytope.vertices((1,1),(2,2),bits=1) == [[1 1]',[0 1]',[1 0]',[0 0]']

        @test LocalPolytope.vertices((2,2),(2,1),bits=1) == [
            [1 1 1 1]',[0 0 1 1]',[1 1 0 0]',[0 0 0 0]'
        ]
        @test issetequal(LocalPolytope.vertices((2,2),(1,2),bits=1), [
            [0 0 0 0 0 0]',[0 0 0 0 1 1]',[0 0 0 0 1 0]',[0 0 0 0 0 1]',
            [1 1 0 0 0 0]',[1 1 1 1 0 0]',[1 1 1 0 0 0]',[1 1 0 1 0 0]',
            [0 1 0 0 0 0]',[0 1 0 1 1 0]',[0 1 0 1 0 0]',[0 1 0 0 1 0]',
            [1 0 0 0 0 0]',[1 0 1 0 0 1]',[1 0 1 0 0 0]',[1 0 0 0 0 1]'
        ])
        @test issetequal(LocalPolytope.vertices((2,1),(2,2),bits=1), [
                [0 0 0 0]',[0 0 0 1]',[0 0 1 0]',[0 0 1 1]',
                [0 1 0 0]',[0 1 0 1]',[0 1 1 0]',[0 1 1 1]',
                [1 0 0 0]',[1 0 0 1]',[1 0 1 0]',[1 0 1 1]',
                [1 1 0 0]',[1 1 0 1]',[1 1 1 0]',[1 1 1 1]'
        ])

        @test issetequal(LocalPolytope.vertices((3,1),(1,2),bits=1), [
            [1 1 1]',[0 0 0]',[1 0 0]',[1 1 0]',[0 1 0]',[0 1 1]',[0 0 1]',[1 0 1]'
        ])
        # there are no strategies which allow for 3 distinct outputs simultaneously
        @test issetequal(LocalPolytope.vertices((3,1),(1,3),bits=1), [
            [1 0 0 0 0 0]',[1 1 0 0 0 0]',[1 0 1 0 0 0]',[1 1 1 0 0 0]',[1 1 0 0 0 1]',[1 0 1 0 1 0]',[1 0 0 0 1 1]',
            [0 0 0 1 0 0]',[0 0 0 1 1 0]',[0 0 0 1 0 1]',[0 1 1 1 0 0]',[0 1 0 1 0 1]',[0 0 1 1 1 0]',[0 0 0 1 1 1]',
            [0 0 0 0 0 0]',[0 1 0 0 0 0]',[0 0 1 0 0 0]',[0 1 1 0 0 0]',[0 0 0 0 1 1]',[0 0 0 0 1 0]',[0 0 0 0 0 1]'
        ])
    end

    @testset "trivial one-way signaling cases" begin
        @test LocalPolytope.vertices((1,1),(1,1),bits=1,fixed_direction=false) == [fill(1,(1,1))]

        @test LocalPolytope.vertices((2,1),(1,1),bits=1,fixed_direction=false) == [[1 1]']
        @test LocalPolytope.vertices((1,2),(1,1),bits=1,fixed_direction=false) == [fill(1,(1,1)),fill(0,(1,1))]
        @test LocalPolytope.vertices((1,1),(2,1),bits=1,fixed_direction=false) == [[1 1]']
        @test LocalPolytope.vertices((1,1),(1,2),bits=1,fixed_direction=false) == [fill(1,(1,1)),fill(0,(1,1))]

        @test LocalPolytope.vertices((2,2),(1,1),bits=1,fixed_direction=false) == [[1 1]',[0 1]',[1 0]',[0 0]']
        @test LocalPolytope.vertices((2,1),(2,1),bits=1,fixed_direction=false) == [[1 1 1 1]']
        @test LocalPolytope.vertices((2,1),(1,2),bits=1,fixed_direction=false) == [[1 1]',[0 0]',[0 1]',[1 0]']
        @test LocalPolytope.vertices((1,2),(2,1),bits=1,fixed_direction=false) == [[1 1]',[0 0]',[0 1]',[1 0]']
        @test LocalPolytope.vertices((1,2),(1,2),bits=1,fixed_direction=false) == [[1 0 0]',[0 1 0]',[0 0 1]',[0 0 0]']
        @test LocalPolytope.vertices((1,1),(2,2),bits=1,fixed_direction=false) == [[1 1]',[0 1]',[1 0]',[0 0]']

        @test issetequal(LocalPolytope.vertices((2,2),(2,1),bits=1,fixed_direction=false), [
            [0 0 0 0]',[0 0 0 1]',[0 0 1 0]',[0 0 1 1]',
            [0 1 0 0]',[0 1 0 1]',[0 1 1 0]',[0 1 1 1]',
            [1 0 0 0]',[1 0 0 1]',[1 0 1 0]',[1 0 1 1]',
            [1 1 0 0]',[1 1 0 1]',[1 1 1 0]',[1 1 1 1]'
        ])
        @test issetequal(LocalPolytope.vertices((2,2),(1,2),bits=1,fixed_direction=false), [
            [1 1 0 0 0 0]',[1 0 0 1 0 0]',[1 0 0 0 0 1]',[1 0 0 0 0 0]',
            [0 1 1 0 0 0]',[0 0 1 1 0 0]',[0 0 1 0 0 1]',[0 0 1 0 0 0]',
            [0 1 0 0 1 0]',[0 0 0 1 1 0]',[0 0 0 0 1 1]',[0 0 0 0 1 0]',
            [0 1 0 0 0 0]',[0 0 0 1 0 0]',[0 0 0 0 0 1]',[0 0 0 0 0 0]'
        ])
        @test issetequal(LocalPolytope.vertices((2,1),(2,2),bits=1,fixed_direction=false), [
            [0 0 0 0]',[0 0 0 1]',[0 0 1 0]',[0 0 1 1]',
            [0 1 0 0]',[0 1 0 1]',[0 1 1 0]',[0 1 1 1]',
            [1 0 0 0]',[1 0 0 1]',[1 0 1 0]',[1 0 1 1]',
            [1 1 0 0]',[1 1 0 1]',[1 1 1 0]',[1 1 1 1]'
        ])
        @test issetequal(LocalPolytope.vertices((1,2),(2,2),bits=1,fixed_direction=false), [
            [1 1 0 0 0 0]',[1 0 0 1 0 0]',[1 0 0 0 0 1]',[1 0 0 0 0 0]',
            [0 1 1 0 0 0]',[0 0 1 1 0 0]',[0 0 1 0 0 1]',[0 0 1 0 0 0]',
            [0 1 0 0 1 0]',[0 0 0 1 1 0]',[0 0 0 0 1 1]',[0 0 0 0 1 0]',
            [0 1 0 0 0 0]',[0 0 0 1 0 0]',[0 0 0 0 0 1]',[0 0 0 0 0 0]'
        ])
    end

    @testset "verify (22)(22) no-signaling setup" begin
        vertices = LocalPolytope.vertices((2,2), (2,2))
        @test size(vertices) == (16,)

        @testset "(22)(22) vertices are well-formed and unique" begin
            for vertex in vertices
                @test size(vertex) == (8,1)
            end

            @test size(vertices) == size(unique(vertices))
        end

        # check results against specialized (22)(22) function
        vertices_match_p8 = map(
            x -> begin
                catv = cat([1],x,dims=1)
                v8 = (Behavior.gen_to_ns_proj((2,2),(2,2))*catv)[2:end]
                reshape(v8,(8,1))
            end,
            DichotomicLocalPolytope.deterministic_no_signaling_behaviors()
        )

        @test issetequal(vertices, vertices_match_p8)
    end

    @testset "verify (22)-1->(22) fixed direction signaling case" begin
        alice_expt = (2,2)
        bob_expt = (2,2)

        vertices = LocalPolytope.vertices(alice_expt, bob_expt, bits=1)
        @test size(vertices) == (64,)
        @test size(vertices) == size(unique(vertices))

        # check against existing results
        vertices_p10 = DichotomicLocalPolytope.deterministic_fixed_direction_signaling_behaviors(subspace="p10")

        # need to recast vertices as matrix and not array
        adjoint_match = []
        for vertex in vertices_p10
            push!(adjoint_match, reshape(vertex, length(vertex),1))
        end
        adjoint_match

        @test issetequal(vertices, adjoint_match)
    end

    @testset "verify (22)<-1->(22) one-way signaling case" begin
        vertices = LocalPolytope.vertices((2,2), (2,2), bits=1, fixed_direction=false)

        @test size(vertices) == (112,)
        @test size(vertices) == size(unique(vertices))

        match_vertices = map(
            (x) -> reshape(x,(12,1)),
            DichotomicLocalPolytope.deterministic_one_way_signaling_behaviors(subspace="p12")
        )

        @test issetequal(vertices, match_vertices)
    end

    test_cases = [
        ((2,2),(2,2),1),
        ((3,2),(3,2),1),
    ]
    @testset "equivalence of dits and bits" for case in test_cases

        vbits = LocalPolytope.vertices(case[1], case[2], bits=case[3])
        vdits = LocalPolytope.vertices(case[1], case[2], dits=2^case[3])

        @test vbits == vdits
    end

    @testset "one dit is the same as no-signaling" for case in test_cases
        vdits = LocalPolytope.vertices(case[1],case[2],dits=1)
        vns = LocalPolytope.vertices(case[1],case[2],rep="fixed-direction")

        @test vdits == vns
    end

    @testset "odd numbered dits" begin
        vdits = LocalPolytope.vertices((3,1),(1,3),dits=3)
        v_match = [
            [0 0 0 0 0 0]',[0 0 1 0 0 0]',[0 1 0 0 0 0]',[0 1 1 0 0 0]',[1 0 0 0 0 0]',
            [1 0 1 0 0 0]',[1 1 0 0 0 0]',[1 1 1 0 0 0]',[0 0 0 0 0 1]',[0 0 0 0 1 0]',
            [0 0 0 0 1 1]',[0 0 0 1 0 0]',[0 0 0 1 0 1]',[0 0 0 1 1 0]',[0 0 0 1 1 1]',
            [0 0 1 0 1 0]',[0 0 1 1 0 0]',[0 0 1 1 1 0]',[0 1 0 0 0 1]',[0 1 0 1 0 0]',
            [0 1 0 1 0 1]',[1 0 0 0 0 1]',[1 0 0 0 1 0]',[1 0 0 0 1 1]',[0 1 1 1 0 0]',
            [1 0 1 0 1 0]',[1 1 0 0 0 1]'
        ]

        @test issetequal(vdits, v_match)
    end
end

@testset "compute_no_signaling_behavior()" begin
    @testset "trivial cases" begin
        @test LocalPolytope.compute_no_signaling_behavior(
            ([1], [fill(1,(1,1))], 1, 1),
            ([1], [fill(1,(1,1))], 1, 1),
            subspace="generalized"
        ) == fill(1,(1,1))
        @test LocalPolytope.compute_no_signaling_behavior(
            ([1 0], [[1 0]',[0 1]'], 2, 1),
            ([1], [fill(1,(1,1))], 1, 1),
            subspace="generalized"
        ) == [1 0]'
    end
end

@testset "strategy_to_behavior()" begin
    @test LocalPolytope.strategy_to_behavior([1 2 3;4 5 6], constant = false) == [1 2 3 4 5 6]'
    @test LocalPolytope.strategy_to_behavior([1 2;3 4;5 6], constant = false) == [1 2 3 4 5 6]'
    @test LocalPolytope.strategy_to_behavior([1 2 3;4 5 6;7 8 9], constant = false) == [1 2 3 4 5 6 7 8 9]'

    @test LocalPolytope.strategy_to_behavior([1 2 3;4 5 6]) == [1 1 2 3 4 5 6]'
    @test LocalPolytope.strategy_to_behavior([1 2;3 4;5 6]) == [1 1 2 3 4 5 6]'
    @test LocalPolytope.strategy_to_behavior([1 2 3;4 5 6;7 8 9]) == [1 1 2 3 4 5 6 7 8 9]'
end

@testset "behavior_to_strategy()" begin
    @test LocalPolytope.behavior_to_strategy(3, 3, [1 1 0 0 0 1 0 0 0 1]') == [1 0 0;0 1 0;0 0 1]
    @test LocalPolytope.behavior_to_strategy(3, 3, [1 0 0 0 1 0 0 0 1]') == [1 0 0;0 1 0;0 0 1]
    @test LocalPolytope.behavior_to_strategy(3, 3, [1 0 0 0 1 1 0 0 0 1]') == [0 0 0;1 1 0;0 0 1]
    @test LocalPolytope.behavior_to_strategy(4, 2, [1 0 1 0 0 1 0 1]') == [1 0 1 0;0 1 0 1]

    @test_throws DomainError LocalPolytope.behavior_to_strategy(3, 3, [1 1 0 0 0 1 0]')
end

@testset "num_prepare_and_measure_vertices()" begin
    for N in 1:5
        for d in 1:5
            α = (N,1)
            β = (1,N)
            vertices = unique(LocalPolytope.vertices(α, β, dits=d))
            @test length(vertices) == LocalPolytope.num_prepare_and_measure_vertices(N,d)
        end
    end
end

@testset "num_inhomogeneous_prepare_and_measure_vertices()" begin
    for X in 1:5
        for B in 1:5
            for d in 1:5
                α = (X,1)
                β = (1,B)
                vertices = unique(LocalPolytope.vertices(α, β, dits=d))
                @test length(vertices) == LocalPolytope.num_inhomogeneous_prepare_and_measure_vertices(X,B,d)
            end
        end
    end
end

@testset "facet_to_matrix()" begin
    @test LocalPolytope.facet_to_matrix(3,3, [-2 1 0 0 0 1 0 0 0 1]) == (2, [1 0 0;0 1 0;0 0 1])
    @test LocalPolytope.facet_to_matrix(3,3, [-1 1 0 -1 0 1 -1 0 0 0]) == (2, [1 0 0;0 1 0;0 0 1])

    @test LocalPolytope.facet_to_matrix(5,5,
        [-1 0 0 -1 1 -1 0 0 0 0 0 0 0 1 -1 -1 1 0 -1 -1 0 0 0 0 0 0]
    ) == (4, [0 0 0 2 0;0 0 1 1 1;0 0 2 0 0;1 0 0 0 1;0 0 1 1 1])
end

end
