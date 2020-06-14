using Test, LinearAlgebra

@testset "(22)<-1->(22) bacon & toner facet inequalities" begin

using BellScenario: LocalPolytope, ConvexPolytope, Degeneracy

function bacon_toner_inequalities()
    inequalities = []

    # 16 positivity inequalities
    for i in 1:16
        ineq = zeros(Int64, (1,17))
        ineq[i+1] = -1
        push!(inequalities, ineq)
    end

    # 16 inequalities of second form
    a_set = [[0 1 0 1], [1 0 1 0], [0 1 1 0], [1 0 0 1]]
    b_set = [[0 0 1 1], [1 1 0 0], [0 1 1 0], [1 0 0 1]]
    for a in a_set
        for b in b_set
            ineq = zeros(Int64, (1,17))

            ineq[1] = -2

            p1_index = a[1]*8 + b[1]*4 + 1 + 1
            p2_index = a[2]*8 + b[2]*4 + 2 + 1
            p3_index = a[3]*8 + b[3]*4 + 3 + 1
            p4_index = a[4]*8 + b[4]*4 + 4 + 1

            ineq[p1_index] = 1
            ineq[p2_index] = 1
            ineq[p3_index] = 1
            ineq[p4_index] = 1

            push!(inequalities, ineq)
        end
    end

    # 16 inequalities of 3rd form
    for a in [0,1]
        for b in [0,1]
            for i in [0,1]
                for j in [0,1]
                    ineq = zeros(Int64, (1,17))

                    p1_index = a*8 + 0*4 + i*2 + j*1 + 2
                    p2_index = a*8 + 1*4 + i*2 + j*1 + 2
                    p3_index = 0*8 + b*4 + (1-i)*2 + (1-j)*1 + 2
                    p4_index = 1*8 + b*4 + (1-i)*2 + (1-j)*1 + 2
                    p5_index = a*8 + b*4 + i*2 + (1-j)*1 + 2

                    ineq[p1_index] = -1
                    ineq[p2_index] = -1
                    ineq[p3_index] = -1
                    ineq[p4_index] = -1
                    ineq[p5_index] = 1

                    push!(inequalities, ineq)
                end
            end
        end
    end

    # Transform Bacon and Toner results to normalization subspace
    transform = [
        1 -1 0 0 0 -1 0 0 0 -1 0 0 0;
        1 0 -1 0 0 0 -1 0 0 0 -1 0 0;
        1 0 0 -1 0 0 0 -1 0 0 0 -1 0;
        1 0 0 0 -1 0 0 0 -1 0 0 0 -1;
    ]

    reduced_inequalities = []
    for ineq in inequalities
        reduced_ineq = ineq[:,14:17]*transform + ineq[:,1:13]
        push!(reduced_inequalities, reduced_ineq)
    end

    reduced_inequalities
end

@testset "verifying polytope" begin

    dir = "./test/regression/files/"

    vertices = LocalPolytope.vertices((2,2), (2,2), bits=1, fixed_direction=false)

    # Bacon and Toner claim 112 unique vertices
    @test size(vertices) == (112,)
    @test size(vertices) == size(unique(vertices))

    constraints = ConvexPolytope.facet_constraints(vertices, dir=dir)
    @test constraints[2] == []
    @test constraints[3] == []

    inequalities = constraints[1]
    # Bacon and Toner claim 48 unique facet inequalities
    @test size(inequalities) == (48,)
    @test size(inequalities) == size(unique(inequalities))

    @test issetequal( inequalities, bacon_toner_inequalities())

    canonical_groups = Degeneracy.canonical_facets((2,2),(2,2),inequalities,"normalized")
    @test length(canonical_groups) == 3
    @test issetequal(inequalities, collect(Iterators.flatten(canonical_groups)))
end

end
