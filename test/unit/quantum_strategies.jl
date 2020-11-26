using Test, QBase

@testset "./src/quantum_strategies.jl" begin

using BellScenario

@testset "quantum_strategy()" begin
    q_strat = quantum_strategy(
        Observables.POVM([[1 0 0;0 0 0;0 0 0],[0 0 0;0 1 0;0 0 0],[0 0 0;0 0 0;0 0 1]]),
        States.DensityMatrix.([[1 0 0;0 0 0;0 0 0],[0 0 0;0 1 0;0 0 0],[0 0 0;0 0 0;0 0 1]])
    )

    @test q_strat isa Strategy
    @test q_strat.scenario == BlackBox(3,3)
    @test q_strat == [1 0 0;0 1 0;0 0 1]
end

@testset "quantum_strategy(LocalSignaling)" begin
    scenario = LocalSignaling(3,3,2)
    q_strat = quantum_strategy(
        Observables.mirror_symmetric_qubit_3povm(π/3),
        States.trine_qubits,
        scenario
    )

    @test q_strat isa Strategy
    @test q_strat ≈ [2/3 1/6 1/6;1/6 2/3 1/6;1/6 1/6 2/3]
    @test q_strat.scenario isa LocalSignaling

    @test_throws DomainError quantum_strategy(
        Observables.POVM([[1 0 0;0 0 0;0 0 0],[0 0 0;0 1 0;0 0 0],[0 0 0;0 0 0;0 0 1]]),
        States.DensityMatrix.([[1 0 0;0 0 0;0 0 0],[0 0 0;0 1 0;0 0 0],[0 0 0;0 0 0;0 0 1]]),
        scenario
    )
end

end
