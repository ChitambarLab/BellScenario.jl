using Test, XPORTA

@testset "chsh violation" begin

using BellScenario: QBase, QuantumBehavior

chsh_ineqs = XPORTA.read_ieq("./test/BellComm/integration/files/22-22_no_signaling.ieq").inequalities
bell_states = (
    1/sqrt(2)*[1;0;0;1],
    1/sqrt(2)*[1;0;0;-1],
    1/sqrt(2)*[0;1;1;0],
    1/sqrt(2)*[0;1;-1;0]
)

@testset "maximal violation with bell states" begin

    povm_a0 = [[1 0; 0 0], [0 0; 0 1]] # z measurement
    povm_a1 = [[0.5 0.5;0.5 0.5], [0.5 -0.5; -0.5 0.5]] # x measurement
    povm_b0 = [
        0.5*([1 0; 0 1] + 1/sqrt(2)*([1 0 ; 0 -1] + [0 1; 1 0])),
        0.5*([1 0; 0 1] - 1/sqrt(2)*([1 0 ; 0 -1] + [0 1; 1 0]))
    ] # -z - x measurement
    povm_b1 = [
        0.5*([1 0; 0 1] + 1/sqrt(2)*([1 0 ; 0 -1] - [0 1; 1 0])),
        0.5*([1 0; 0 1] - 1/sqrt(2)*([1 0 ; 0 -1] - [0 1; 1 0]))
    ] # z - x measuremnt

    @test QBase.is_povm(povm_a0)
    @test QBase.is_povm(povm_a1)
    @test QBase.is_povm(povm_b0)
    @test QBase.is_povm(povm_b1)

    for ψ in bell_states
        behavior = QuantumBehavior.bipartite_scenario([povm_a0,povm_a1],[povm_b0,povm_b1],ψ)

        violation_count = 0
        for ineq in eachrow(chsh_ineqs)

            exp_val = (cat(-1*ineq[end],ineq[1:end-1]',dims=2)*behavior)[1]
            if (exp_val > 0) & !(exp_val + 1 ≈ 1)
                violation_count += 1
            end
        end

        @test violation_count == 1
    end
end
end
