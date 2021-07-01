using Test, QBase

@testset "./src/quantum_strategies.jl" begin

using BellScenario

@testset "quantum_strategy()" begin
    @testset "QBase types" begin
        q_strat = quantum_strategy(
            POVM([[1 0 0;0 0 0;0 0 0],[0 0 0;0 1 0;0 0 0],[0 0 0;0 0 0;0 0 1]]),
            State.([[1 0 0;0 0 0;0 0 0],[0 0 0;0 1 0;0 0 0],[0 0 0;0 0 0;0 0 1]])
        )

        @test q_strat isa Strategy
        @test q_strat.scenario == BlackBox(3,3)
        @test q_strat == [1 0 0;0 1 0;0 0 1]
    end

    @testset "matrix types" begin
        q_strat = quantum_strategy(
            [[1 0;0 0],[0 0;0 1]],
            [[1 0;0 0],[0 0;0 1]]
        )

        @test q_strat isa Strategy
        @test q_strat.scenario == BlackBox(2,2)
        @test q_strat == [1 0;0 1]
    end
end

@testset "quantum_strategy(LocalSignaling)" begin
    @testset "QBase types" begin
        scenario = LocalSignaling(3,3,2)
        q_strat = quantum_strategy(
            mirror_symmetric_qubit_3povm(π/3),
            trine_qubit_states(),
            scenario
        )

        @test q_strat isa Strategy
        @test q_strat ≈ [2/3 1/6 1/6;1/6 2/3 1/6;1/6 1/6 2/3]
        @test q_strat.scenario isa LocalSignaling
    end

    @testset "Matrix types" begin
        scenario = LocalSignaling(2,2,2)

        q_strat = quantum_strategy(
            [[1 0;0 0],[0 0;0 1]],
            [[1 0;0 0],[0 0;0 1]],
            scenario
        )

        @test q_strat isa Strategy
        @test q_strat.scenario == LocalSignaling(2,2,2)
        @test q_strat == [1 0;0 1]
    end

    @testset "DomainErrors" begin
        scenario = LocalSignaling(3,3,2)
        @test_throws DomainError quantum_strategy(
            POVM([[1 0 0;0 0 0;0 0 0],[0 0 0;0 1 0;0 0 0],[0 0 0;0 0 0;0 0 1]]),
            State.([[1 0 0;0 0 0;0 0 0],[0 0 0;0 1 0;0 0 0],[0 0 0;0 0 0;0 0 1]]),
            scenario
        )
    end
end

@testset "quantum_strategy(BipartiteNonSignaling)" begin
    @testset "QBase types ch violation" begin
        scenario = BipartiteNonSignaling(2,2,2,2)
        ρ_AB = State([1 0 0 1;0 0 0 0;0 0 0 0;1 0 0 1]/2)
        Π_Ax = [
            POVM([[1. 0;0 0],[0 0;0 1.]]),
            POVM([[1 1;1 1]/2,[1 -1;-1 1]/2])
        ]
        Π_By = [
            POVM([bloch_qubit_state(π/4),bloch_qubit_state(-3π/4)]),
            POVM([bloch_qubit_state(3π/4),bloch_qubit_state(-π/4)])
        ]

        q_strat = quantum_strategy(ρ_AB, Π_Ax, Π_By, scenario)
        @test q_strat isa Strategy
        @test q_strat.scenario isa BipartiteNonSignaling

        # generalized rep of CH inequality
        g = BellGame([1 0 1 1;1 1 0 0;0 1 1 0;1 1 1 0], 3)
        @test sum(g .* q_strat) ≈ 3.207 atol=2e-4
    end

    @testset "matrix types classical states" begin
        scenario = BipartiteNonSignaling(2,2,2,2)
        ρ_AB = kron([1 0;0 0],[0 0;0 1])
        Π_Ax = [
            [[1 0;0 0],[0 0;0 1]],
            [[0 0;0 1],[1 0;0 0]]
        ]
        Π_By = [
            [[1 0;0 0],[0 0;0 1]],
            [[0 0;0 1],[1 0;0 0]]
        ]

        q_strat = quantum_strategy(ρ_AB, Π_Ax, Π_By, scenario)
        @test q_strat isa Strategy
        @test q_strat == [0 1 0 0;1 0 0 0;0 0 0 1;0 0 1 0]
        @test q_strat.scenario == scenario
    end

    @testset "Domain Errors" begin
        ρ_AB = State(kron([1 0;0 0],[0 0;0 1]))
        Π = POVM([[1 0;0 0],[0 0;0 1]])
        Π_Ax = [Π, Π]
        Π_By = [Π, Π]

        scenario = BipartiteNonSignaling(3,2,2,2)
        @test_throws DomainError quantum_strategy(ρ_AB, Π_Ax, Π_By, scenario)

        scenario = BipartiteNonSignaling(2,3,2,2)
        @test_throws DomainError quantum_strategy(ρ_AB, Π_Ax, Π_By, scenario)

        scenario = BipartiteNonSignaling(2,2,3,2)
        @test_throws DomainError quantum_strategy(ρ_AB, Π_Ax, Π_By, scenario)

        scenario = BipartiteNonSignaling(2,2,2,3)
        @test_throws DomainError quantum_strategy(ρ_AB, Π_Ax, Π_By, scenario)
    end
end

end
