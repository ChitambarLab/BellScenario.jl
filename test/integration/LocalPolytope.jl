using Test

@testset "/src/LocalPolytope.jl" begin

using BellScenario: LocalPolytope

@testset "permutation encodings are not necessary for communication protocols" begin

    # this test verifies that it is okay to use the communication protocols method
    # instead of the strategies method to compute communication protocols.

    # testing for all 1-bit cases with 3 in/out or fewer
    @testset "bob input $y" for y in 1:3
        @testset "bob output $b" for b in 1:3
            @testset "alice input = $x" for x in 1:3
                @testset "protocol out = $num_out" for num_out in 2:2
                    I = QMath.id(y)

                    ρ_perms = map(
                        protocol -> kron(I,protocol),
                        LocalPolytope.strategies(x,num_out)
                    )
                    ρ_test = map(
                        protocol -> kron(I,protocol),
                        LocalPolytope.communication_protocols(x,num_out)
                    )

                    β_joint_set = LocalPolytope.strategies(num_out*y,b)


                    test_set = []
                    perm_set = []
                    for ρ in ρ_test
                        test_set = cat(test_set, map(
                            (β) -> β*ρ,
                            β_joint_set
                        ), dims = 1)
                    end

                    for ρ in ρ_perms
                        perm_set = cat(perm_set, map(
                            (β) -> β*ρ,
                            β_joint_set
                        ), dims = 1)
                    end

                    @test issetequal(unique(test_set), unique(perm_set))
                end
            end
        end
    end
end

end
