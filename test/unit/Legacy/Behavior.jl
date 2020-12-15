using Test, LinearAlgebra

@testset "/src/Legacy/Behavior.jl" begin

using BellScenario: Behavior, QMath

@testset "gen_to_fd_proj()" begin
    @testset "dichotomic case" begin
        @test Behavior.gen_to_fd_proj((2,2),(2,2)) == [
            1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0;
            0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
        ]
    end
end

@testset "fd_to_gen_proj()" begin
    @testset "dichotomic case" begin
        @test Behavior.fd_to_gen_proj((2,2),(2,2)) == [
            1 0 0 0 0 0 0 0 0 0 0;
            0 0 0 1 0 0 0 0 0 0 0;
            0 0 0 0 1 0 0 0 0 0 0;
            0 0 0 0 0 1 0 0 0 0 0;
            0 0 0 0 0 0 1 0 0 0 0;
            0 1 0 -1 0 0 0 0 0 0 0;
            0 1 0 0 -1 0 0 0 0 0 0;
            0 0 1 0 0 -1 0 0 0 0 0;
            0 0 1 0 0 0 -1 0 0 0 0;
            0 0 0 0 0 0 0 1 0 0 0;
            0 0 0 0 0 0 0 0 1 0 0;
            0 0 0 0 0 0 0 0 0 1 0;
            0 0 0 0 0 0 0 0 0 0 1;
            1 -1 0 0 0 0 0 -1 0 0 0;
            1 -1 0 0 0 0 0 0 -1 0 0;
            1 0 -1 0 0 0 0 0 0 -1 0;
            1 0 -1 0 0 0 0 0 0 0 -1;
        ]
    end

end

@testset "ns_to_gen_proj()" begin
    @testset "dichotomic case" begin
        @test Behavior.ns_to_gen_proj((2,2),(2,2)) == [
            1 0 0 0 0 0 0 0 0;
            0 0 0 0 0 1 0 0 0;
            0 0 0 0 0 0 1 0 0;
            0 0 0 0 0 0 0 1 0;
            0 0 0 0 0 0 0 0 1;
            0 1 0 0 0 -1 0 0 0;
            0 1 0 0 0 0 -1 0 0;
            0 0 1 0 0 0 0 -1 0;
            0 0 1 0 0 0 0 0 -1;
            0 0 0 1 0 -1 0 0 0;
            0 0 0 0 1 0 -1 0 0;
            0 0 0 1 0 0 0 -1 0;
            0 0 0 0 1 0 0 0 -1;
            1 -1 0 -1 0 1 0 0 0;
            1 -1 0 0 -1 0 1 0 0;
            1 0 -1 -1 0 0 0 1 0;
            1 0 -1 0 -1 0 0 0 1;
        ]
    end

    test_cases = [
        ((2,2),(2,2),(17,9)),
        ((2,3),(2,3),(37,25)),
        ((3,2),(2,2),(25,12)),
        ((3,2),(3,2),(37,16)),
        ((3,3),(3,3),(82,49)),
        ((4,2),(4,2),(65,25)),
        ((2,3),(4,4),(97,65)),
        ((4,4),(4,4),(257,169))
    ]

    @testset "dimensionality $(case[1]), $(case[2])" for case in test_cases
        @test size(Behavior.ns_to_gen_proj(case[1],case[2])) == case[3]
    end

    @testset "behavior validity $(case[1]), $(case[2])" for case in test_cases

        α_expt = case[1]
        β_expt = case[2]
        α_dim = α_expt[1]*(α_expt[2] - 1)
        β_dim = β_expt[1]*(β_expt[2] - 1)

        α_val = round(1/α_expt[2],sigdigits = 3)
        β_val = round(1/β_expt[2],sigdigits = 3)
        αβ_val = round(α_val*β_val, sigdigits = 3)

        ns_p = cat([1],fill(α_val,α_dim),fill(β_val,β_dim), fill(αβ_val,α_dim*β_dim),dims =1)
        @test Behavior.is_valid(ns_p, α_expt, β_expt, "no-signaling")

        gen_p = Behavior.ns_to_gen_proj(α_expt,β_expt)*ns_p
        @test Behavior.is_valid(gen_p,case[1],case[2])
    end

    @testset "invertibility $(case[1]), $(case[2])" for case in test_cases
        ns = Behavior.gen_to_ns_proj(case[1],case[2])
        ns_inv = Behavior.ns_to_gen_proj(case[1],case[2])

        @test diagm( 0=> fill(1,(case[3])[2])) == ns*ns_inv
    end
end

@testset "norm_to_gen_proj()" begin
    @testset "dichotomic setup" begin
        proj_match = cat(
            diagm( 0 => fill(1, 13) ),
            [
                1 -1 0 0 0 -1 0 0 0 -1 0 0 0;
                1 0 -1 0 0 0 -1 0 0 0 -1 0 0;
                1 0 0 -1 0 0 0 -1 0 0 0 -1 0;
                1 0 0 0 -1 0 0 0 -1 0 0 0 -1;
            ],
            dims = 1
        )

        Behavior.norm_to_gen_proj((2,2),(2,2))
        @test proj_match == Behavior.norm_to_gen_proj((2,2),(2,2))
    end
end

@testset "gen_to_norm_proj()" begin
    @testset "dichotomic scenario" begin
        proj_match = cat(diagm( 0 => fill(1,13)),zeros((13,4)),dims=2)
        @test proj_match == Behavior.gen_to_norm_proj((2,2),(2,2))
    end

    @testset "3333 scenario" begin
        proj_match = cat(diagm( 0 => fill(1,73)),zeros(73,9),dims=2)
        @test proj_match == Behavior.gen_to_norm_proj((3,3),(3,3))
    end
end

@testset "gen_to_ns_proj()" begin
    test_cases = [
        (
            [1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;0;0],
            (2,2),(2,2),
            [1;1;1;1;1;1;1;1;1]
        ),
        (
            [1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;],
            (2,3),(2,3),
            [1,1,1,0,0,1,1,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0]
        ),
        (
            [1;0;0;0;0;0;0;0;0;1;1;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;],
            (2,3),(2,3),
            [1;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0],
        )
    ]

    @testset "maps to valid no-signaling behavior" for case in test_cases
        p = case[1]
        α_expt = case[2]
        β_expt = case[3]
        p_match = case[4]

        @test Behavior.is_valid(p,α_expt,β_expt)
        @test Behavior.is_valid(p_match, α_expt, β_expt, "no-signaling")

        p_ns = Behavior.gen_to_ns_proj(α_expt,β_expt)*p

        @test Behavior.is_valid(p_ns, α_expt, β_expt, "no-signaling")
        @test p_ns ==  p_match
    end
end

@testset "index()" begin
    @testset "rep = generalized" begin
        α_expt = (2,2)
        β_expt = (2,2)

        match_val = 2
        for a in 1:α_expt[2]
            for b in 1:β_expt[2]
                for x in 1:α_expt[1]
                    for y in 1:β_expt[1]
                        id = Behavior.index((a,b,x,y), α_expt, β_expt)
                        @test id == match_val
                        match_val += 1
                    end
                end
            end
        end

        α_expt = (5,2)
        β_expt = (3,4)

        @test Behavior.index((1,1,1,1),α_expt,β_expt) == 2
        @test Behavior.index((2,4,5,3),α_expt,β_expt) == 121
        @test Behavior.index((1,1,5,1),α_expt,β_expt) == 14

        α_expt = (3,3)
        β_expt = (1,3)

        @test Behavior.index((1,1,1,1),α_expt,β_expt) == 2
        @test Behavior.index((3,3,3,1),α_expt,β_expt) == 28
        @test Behavior.index((1,1,3,1),α_expt,β_expt) == 4
        @test Behavior.index((3,1,1,1),α_expt,β_expt) == 20
        @test Behavior.index((1,3,1,1),α_expt,β_expt) == 8
    end

    @testset "rep = normalized" begin
        α_expt = (2,2)
        β_expt = (2,2)

        @test Behavior.index((1,1,1,1),α_expt,β_expt,"normalized") == 2
        @test Behavior.index((1,2,1,1),α_expt,β_expt,"normalized") == 6
        @test Behavior.index((2,1,1,1),α_expt,β_expt,"normalized") == 10

        for x in 1:2
            for y in 1:2
                @test_throws ArgumentError Behavior.index((2,2,x,y),α_expt,β_expt,"normalized") == 13
            end
        end
    end

    @testset "rep = no-signaling" begin
        α_expt = (2,2)
        β_expt = (2,2)

        @test Behavior.index((1,1,1,1),α_expt,β_expt,"no-signaling") == 6
        @test Behavior.index((0,1,0,1),α_expt,β_expt,"no-signaling") == 4
        @test Behavior.index((1,0,1,0),α_expt,β_expt,"no-signaling") == 2
        @test Behavior.index((2,2,1,1),(2,3),(2,3),"no-signaling") == 22


        @test_throws ArgumentError Behavior.index((2,1,1,0),α_expt,β_expt,"no-signaling")
        @test_throws ArgumentError Behavior.index((1,0,0,1),α_expt,β_expt,"no-signaling")
        @test_throws ArgumentError Behavior.index((0,1,1,0),α_expt,β_expt,"no-signaling")
    end

    @testset "rep == fixed-direction" begin
        @test Behavior.index((1,1,1,1),(2,2),(2,2),"fixed-direction") == 4
        @test Behavior.index((1,0,1,0),(2,2),(2,2),"fixed-direction") == 2
        @test Behavior.index((2,1,2,2),(2,2),(2,2),"fixed-direction") == 11
        @test Behavior.index((2,1,3,2),(3,2),(2,2),"fixed-direction") == 16

        @test_throws ArgumentError Behavior.index((1,2,2,2),(2,2),(2,2),"fixed-direction")
        @test_throws ArgumentError Behavior.index((0,1,0,1),(2,2),(2,2),"fixed-direction")
    end
end

@testset "is_valid()" begin

    gen_cases =[
        ([1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;0;0],(2,2),(2,2), true),
        ([1;.25;.25;.25;.25;.25;.25;.25;.25;.25;.25;.25;.25;.25;.25;.25;.25],(2,2),(2,2), true),
    ]

    @testset "generalized" for case in gen_cases
        @test Behavior.is_valid(case[1],case[2],case[3]) == case[4]
    end

    ns_cases = [
        ([1;1;1;1;1;1;1;1;1],(2,2),(2,2), true),
        ([1;1;0;0;1;0;1;0;0],(2,2),(2,2), true),
        ([1;1;1;0;0;1;1;0;0], (2,2),(2,2), false),
        ([1;0;1;1;0;0;0;1;0], (2,2),(2,2), true),
        ([1;1;1;0;0;0;0;1;1;0;0;0;0;1;1;1;1;0;0;0;0;0;0;0;0], (2,3),(2,3),true),

    ]

    @testset "no-signaling" for case in ns_cases
        @test Behavior.is_valid(case[1],case[2],case[3],"no-signaling") == case[4]
    end

    fd_cases = [
        ([1;1;1;1;1;0;0],(2,2),(1,2),true),
        ([1;1;1;1;0;0;0],(2,2),(1,2),true),
        ([1;0;0;0;0;0;0;0;1;1;0],(2,2),(2,2),true),
        (1/4*[1;2;2;1;1;1;1;1;1;1;1],(2,2),(2,2),true),
        ([1;0;0;0;1;1;0;0;0;0;0;0;0;1;1;0;0;0;0;0;0;1;0;0;0],(3,3),(1,3),true),
        ([1;1;0;0;1;0;0;0;0;0;0;0;0;0;0;0],(3,2),(2,2),true)
    ]

    @testset "fixed-direction" for case in fd_cases
        @test Behavior.is_valid(case[1],case[2],case[3],"fixed-direction") == case[4]
    end

end

@testset "dimension()" begin
    test_cases = [
        ((2,2),(2,2),9,11,13,17),
        ((2,3),(2,3),25,29,33,37),
        ((3,2),(2,2),12,16,19,25),
        ((3,2),(3,2),16,22,28,37),
        ((3,3),(3,3),49,61,73,82),
        ((4,2),(4,2),25,37,49,65),
        ((2,3),(4,4),65,77,89,97),
        ((4,4),(4,4),169,205,241,257)
    ]

    @testset "$(case[1]), $(case[2])" for case in test_cases
        @test Behavior.dimension(case[1],case[2],"generalized") == case[6]
        @test Behavior.dimension(case[1],case[2],"normalized") == case[5]
        @test Behavior.dimension(case[1],case[2],"fixed-direction") == case[4]
        @test Behavior.dimension(case[1],case[2],"no-signaling") == case[3]
    end
end

@testset "conditionals()" begin
    conditionals = Behavior.conditionals((3,1),(1,3),[1;0;1;0;1;0],rep="normalized")
    @test conditionals == [1 0 0;0 1 0;1 0 0]'
    @test conditionals isa QMath.Conditionals

    conditionals = Behavior.conditionals((3,1),(1,3),[1;1;0;1;0;1;0],rep="normalized")
    @test conditionals == [1 0 0;0 1 0;1 0 0]'
    @test conditionals isa QMath.Conditionals

    @test_throws ArgumentError Behavior.conditionals((3,1),(1,3),[1;1;0;1;1;1;0],rep="normalized")
end

@testset "projector algebra" begin

    test_cases = [
        ((2,2),(2,2)),
        ((3,3),(1,3)),
        ((3,2),(3,2))
    ]

    @testset "gen <-> norm α = $(case[1]), β = $(case[2])" for case in test_cases
        α = case[1]
        β = case[2]

        num_out = α[2]*β[2]

        norm_dim = Behavior.dimension(α,β,"normalized")
        norm_id = diagm( 0 => fill(1,norm_dim))

        gen_dim = Behavior.dimension(α,β,"generalized")
        gen_behavior = cat(1,fill(1/num_out, (gen_dim-1)),dims=1 )
        @test Behavior.is_valid(gen_behavior, α, β, "generalized")

        proj = Behavior.gen_to_norm_proj(α,β)
        proj_inv = Behavior.norm_to_gen_proj(α,β)

        @test proj*proj_inv == norm_id
        @test (proj_inv*proj)*gen_behavior ≈ gen_behavior
    end

    @testset "gen <-> fd α = $(case[1]), β = $(case[2])" for case in test_cases
        α = case[1]
        β = case[2]

        num_out = α[2]*β[2]

        fd_dim = Behavior.dimension(α,β,"fixed-direction")
        fd_id = diagm( 0 => fill(1,fd_dim))

        gen_dim = Behavior.dimension(α,β,"generalized")
        gen_behavior = cat(1,fill(1/num_out, (gen_dim-1)),dims=1 )

        @test Behavior.is_valid(gen_behavior, α, β, "generalized")

        proj = Behavior.gen_to_fd_proj(α,β)
        proj_inv = Behavior.fd_to_gen_proj(α,β)

        @test proj*proj_inv == fd_id
        @test (proj_inv*proj)*gen_behavior ≈ gen_behavior
    end

end

@testset "has_constant()" begin
    @test !Behavior.has_constant((2,2),(2,2),[0;0;0;0;0;0;0;0;0;0;0;0])
    @test Behavior.has_constant((2,2),(2,2),[1;0;0;0;0;0;0;0;0;0;0;0;0])
    @test Behavior.has_constant((2,2),(2,2),[1;0;0;0;0;0;0;0;0], rep="no-signaling")

    @test_throws ArgumentError Behavior.has_constant((2,2),(2,2),[0;0;0;0;0;0;0;0;0;0;0;0], rep="generalized")
end

@testset "add_constant()" begin
    behavior1 = [1;0;0;0;1;0]
    behavior2 = [1;0;0;0;0;0;0]

    @test Behavior.add_constants((3,1),(1,3), [behavior1, behavior2]) == [[1;1;0;0;0;1;0],[1;0;0;0;0;0;0]]
end

@testset "remove_constant()" begin
    behavior1 = [1;0;0;0;0;0;0]
    behavior2 = [1;0;0;0;1;0]

    @test Behavior.remove_constants((3,1),(1,3), [behavior1, behavior2]) == [[0;0;0;0;0;0],[1;0;0;0;1;0]]
end

@testset "behavior projection invertibility" begin
    test_cases = [
        ((2,2),(2,2)),
        ((2,3),(2,3)),
        ((3,2),(2,2)),
        ((3,2),(3,2)),
        ((3,3),(3,3)),
        ((4,2),(4,2)),
        ((2,3),(4,4)),
        ((4,4),(4,4))
    ]

    @testset "case: $(case[1]), $(case[2])" for case in test_cases
        ns_proj = Behavior.gen_to_ns_proj(case[1],case[2])
        ns_inv = Behavior.ns_to_gen_proj(case[1],case[2])
        ns_dim = Behavior.dimension(case[1],case[2],"no-signaling")
        @test ns_proj*ns_inv == diagm(0 => fill(1,ns_dim))

        fd_proj = Behavior.gen_to_fd_proj(case[1],case[2])
        fd_inv = Behavior.fd_to_gen_proj(case[1],case[2])
        fd_dim = Behavior.dimension(case[1],case[2],"fixed-direction")

        @test fd_proj*fd_inv == diagm(0 => fill(1,fd_dim))
    end
end

end
