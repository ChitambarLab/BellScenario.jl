using Test

@testset "CHSH facet regression using relabels" begin

using BellScenario: Degeneracy, LocalPolytope, ConvexPolytope

function chsh_inequalities()
    # positivity
    I_ineq = [0 0 0 0 0 -1 0 0 0]
    # non-trivial CH inequality
    II_ineq = [0 -1 0 -1 0 1 1 1 -1]

    (α_out_relabels, β_out_relabels) = Degeneracy.bipartite_output_relabels((2,2),(2,2),"no-signaling")
    in_relabels = Degeneracy.bipartite_input_relabels((2,2),(2,2),"no-signaling")

    ineqs = []
    # TODO: resolve the 3x redundancy in relabeled inequalities
    for α_out_relabel in α_out_relabels
        for β_out_relabel in β_out_relabels
            for in_relabel in in_relabels
                push!(ineqs, I_ineq*in_relabel*α_out_relabel*β_out_relabel)
                push!(ineqs, II_ineq*in_relabel*α_out_relabel*β_out_relabel)
            end
        end
    end

    unique(ineqs)
end

@testset "verifying polytope" begin

    dir = "./test/regression/files/"

    vertices = LocalPolytope.vertices((2,2),(2,2))
    @test size(vertices) == (16,)
    @test size(vertices[1]) == (8,1)

    constraints = ConvexPolytope.facet_constraints(vertices, dir=dir)
    @test constraints[2] == []
    @test constraints[3] == []

    ineq_matches = chsh_inequalities()

    inequalities = constraints[1]
    issetequal(inequalities, ineq_matches)
end

end
