using Test, LinearAlgebra

@testset "test/integration/facet_degeneracies.jl" begin

using BellScenario: Degeneracy

@testset "31->13 polytope" begin
    α = (3,1)
    β = (1,3)

    facets = [
        [0 -1 0 0 0 0 0], [0 0 -1 0 0 0 0], [0 0 0 -1 0 0 0],
        [0 0 0 0 -1 0 0], [0 0 0 0 0 -1 0], [0 0 0 0 0 0 -1],
        [-1 0 0 1 0 0 1], [-1 0 1 0 0 1 0], [-1 1 0 0 1 0 0],
        [-1 -1 0 1 -1 1 0], [-1 -1 1 0 -1 0 1], [-1 0 -1 1 1 -1 0],
        [-1 0 1 -1 1 0 -1], [-1 1 -1 0 0 -1 1], [-1 1 0 -1 0 1 -1]
    ]

    out_sym1 = facets[[1,4,9]]
    out_sym2 = facets[[2,5,8]]
    out_sym3 = facets[[3,6,7]]
    out_sym4 = facets[10:15]

    in_sym1 = facets[1:9]
    in_sym2 = facets[10:15]

    input_relabels = Degeneracy.bipartite_input_relabels(α,β,"normalized")
    output_relabels = Degeneracy.bipartite_output_relabels(α,β,"normalized")

    @testset "verify output relabeling symmetry groups" for sym_group in [out_sym1, out_sym2, out_sym3, out_sym4]
        for i in [1,2]
            for relabel in output_relabels[i]
                for facet in sym_group
                    @test in(facet*relabel,sym_group)
                end
            end
        end
    end

    @testset "verify symmetry groups 1,2,3 are closed under input relabelings" for sym_group in [in_sym1, in_sym2]
        for relabel in input_relabels
            for facet in sym_group
                @test in(facet*relabel,sym_group)
            end
        end
    end
end

end
