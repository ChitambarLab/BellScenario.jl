using Test

@testset "(32)-1->(22) fixed-direction signaling regression" begin

using BellScenario: Degeneracy, LocalPolytope, ConvexPolytope

# """
# maxwell_chitambar_inequalities()
#     Uses input/output relabelings to construct the complete set of bell inequalities
#     from a canonical subset of base inequalities.
#
# Outputs:
#     ineqs: Array, contains inequality vectors.
# """
function maxwell_chitambar_inequalities()

    joint_pos_ineq = [0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0]
    marginal_pos_ineq = [0 -1 0 0 1 0 0 0 0 0 0 0 0 0 0 0]
    I_ineq = [-1 -1 0 0 0 0 -1 1 1 0 -1 -1 0 0 1 0]
    II_ineq = [-1 0 0 0 -1 -1 -1 1 1 0 -1 -1 0 0 1 0]
    III_ineq = [-1 -1 -1 0 0 0 1 1 -1 1 -1 -1 0 0 1 -1]
    IV_ineq = [-1 0 0 0 -1 -1 -1 1 1 0 -1 -1 -1 1 1 0]
    V_ineq = [-1 -1 0 0 1 1 -1 -1 -1 1 0 0 -1 -1 1 -1]
    VI_ineq = [-1 -3 0 0 2 2 -2 1 0 -1 -1 -1 -1 0 1 -2]
    VII_ineq = [-1 -3 0 0 2 2 -2 1 0 -1 -1 -1 -2 1 1 -2]
    VIII_ineq = [-1 -3 0 0 2 2 -2 1 1 -2 -1 -1 -2 1 1 -2]

    ineqs = []

    # TODO: construct relabelngs in a way that is inclusive, but does not overshoot
    # by a factor of 4.5x duplicates.
    in_relabels = Degeneracy.bipartite_input_relabels((3,2),(2,2),"fixed-direction")
    (α_out_relabels, β_out_relabels) = Degeneracy.bipartite_output_relabels((3,2),(2,2),"fixed-direction")
    for in_relabel in in_relabels
        for α_out_relabel in α_out_relabels
            for β_out_relabel in β_out_relabels

                relabel = in_relabel*α_out_relabel*β_out_relabel

                push!(ineqs, joint_pos_ineq*relabel)
                push!(ineqs, marginal_pos_ineq*relabel)
                push!(ineqs, I_ineq*relabel)
                push!(ineqs, II_ineq*relabel)
                push!(ineqs, III_ineq*relabel)
                push!(ineqs, IV_ineq*relabel)
                push!(ineqs, V_ineq*relabel)
                push!(ineqs, VI_ineq*relabel)
                push!(ineqs, VII_ineq*relabel)
                push!(ineqs, VIII_ineq*relabel)
            end
        end
    end

    # TODO: fix so that we do not need to de-dupe
    unique(ineqs)
end

@testset "verifying polytope" begin

    dir = "./test/regression/files/"

    vertices = LocalPolytope.vertices((3,2),(2,2), bits=1)

    @test size(vertices) == (320,)
    @test size(vertices[1]) == (15,1)

    constraints = ConvexPolytope.facet_constraints(vertices, dir=dir)
    @test constraints[2] == []
    @test constraints[3] == []

    ineq_matches = maxwell_chitambar_inequalities()

    inequalities = constraints[1]
    @test issetequal(inequalities, ineq_matches)
end

end
