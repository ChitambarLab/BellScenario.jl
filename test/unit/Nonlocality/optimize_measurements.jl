using Test, QBase

@testset "./src/Nonlocality/optimize_measurements.jl" begin

using BellScenario

@testset "Nonlocality.optimize_measurement(LocalSignaling)" begin
    @testset "trine states" begin
        scenario = LocalSignaling(3,3,2)
        game = BellGame([1 0 0;0 1 0;0 0 1], 2)
        ρ_states = States.trine_qubits

        dict = Nonlocality.optimize_measurement(scenario, game, ρ_states)

        @test isapprox(dict["violation"], 0.0, atol=1e-6)
        @test all(isapprox.(dict["povm"][1], 2/3*ρ_states[1], atol=1e-6))
        @test all(isapprox.(dict["povm"][2], 2/3*ρ_states[2], atol=1e-6))
        @test all(isapprox.(dict["povm"][3], 2/3*ρ_states[3], atol=1e-6))

        @test dict["states"] == ρ_states
        @test dict["game"] == game
        @test dict["scenario"] == scenario
    end

    @testset "bb84 states" begin
        scenario = LocalSignaling(4,4,2)
        game = BellGame([1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1], 2)
        ρ_states = States.bb84_qubits

        dict = Nonlocality.optimize_measurement(scenario, game, ρ_states)

        @test isapprox(dict["violation"], 0.0, atol=1e-6)
        @test all(isapprox.(dict["povm"][1], 1/2*ρ_states[1], atol=1e-3))
        @test all(isapprox.(dict["povm"][2], 1/2*ρ_states[2], atol=1e-3))
        @test all(isapprox.(dict["povm"][3], 1/2*ρ_states[3], atol=1e-3))
        @test all(isapprox.(dict["povm"][4], 1/2*ρ_states[4], atol=1e-3))
    end

    @testset "sic qubit states" begin
        scenario = LocalSignaling(4,4,2)
        game = BellGame([1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1], 2)
        ρ_states = States.sic_qubits

        dict = Nonlocality.optimize_measurement(scenario, game, ρ_states)

        @test isapprox(dict["violation"], 0.0, atol=1e-5)
        @test all(isapprox.(dict["povm"][1], 1/2*ρ_states[1], atol=1e-5))
        @test all(isapprox.(dict["povm"][2], 1/2*ρ_states[2], atol=1e-5))
        @test all(isapprox.(dict["povm"][3], 1/2*ρ_states[3], atol=1e-5))
        @test all(isapprox.(dict["povm"][4], 1/2*ρ_states[4], atol=1e-5))
    end

    @testset "Errors" begin
        scenario = LocalSignaling(3,3,2)
        states = States.bb84_qubits
        game = BellGame([1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1], 2)

        @test_throws DomainError Nonlocality.optimize_measurement(scenario, game, states)

        scenario = LocalSignaling(4,4,3)
        @test_throws DomainError  Nonlocality.optimize_measurement(scenario, game, states)
    end
end

@testset "Nonlocality.optimize_measurement(::BipartiteNonSignaling)" begin
    @testset "CHSH inequality" begin
        scenario = BipartiteNonSignaling(2,2,2,2)
        game = BellGame([0 1 1 0;0 0 0 1;0 0 0 1;1 0 0 1],2)
        ρ_AB = States.bell_states[1]
        POVMs = [
            Observables.POVM([[1 0;0 0],[0 0;0 1]]),
            Observables.POVM([[0.5 0.5;0.5 0.5],[0.5 -0.5;-0.5 0.5]])
        ]

        @test game.β == 2

        opt_dictA = Nonlocality.optimize_measurement(scenario, game, ρ_AB, A_POVMs=POVMs)
        opt_dictB = Nonlocality.optimize_measurement(scenario, game, ρ_AB, B_POVMs=POVMs)

        @test length(keys(opt_dictA)) == 7
        @test opt_dictA["scenario"] == scenario
        @test opt_dictA["game"] == game
        @test opt_dictA["game"].β == 2
        @test opt_dictA["score"] ≈ 2.20710279
        @test opt_dictA["violation"] ≈ 0.207102796
        @test opt_dictA["state"] == ρ_AB
        @test opt_dictA["A_POVMs"] == POVMs
        @test isapprox(
            opt_dictA["B_POVMs"],
            [
                [States.bloch_qubit(π/4,0), States.bloch_qubit(5π/4,0)],
                [States.bloch_qubit(7π/4,0), States.bloch_qubit(3π/4,0)]
            ],
            atol=1e-5
        )
        @test opt_dictA["B_POVMs"] isa Vector{Observables.POVM}

        @test length(keys(opt_dictB)) == 7
        @test opt_dictB["scenario"] == scenario
        @test opt_dictB["game"] == game
        @test opt_dictB["game"].β == 2
        @test opt_dictB["score"] ≈ 2.2071028
        @test opt_dictB["violation"] ≈ 0.207102796
        @test opt_dictB["state"] == ρ_AB
        @test opt_dictB["B_POVMs"] == POVMs
        @test isapprox(
            opt_dictB["A_POVMs"],
            [
                [States.bloch_qubit(π/4,0), States.bloch_qubit(5π/4,0)],
                [States.bloch_qubit(7π/4,0),States.bloch_qubit(3π/4,0)]
            ],
            atol=1e-5
        )
        @test opt_dictB["A_POVMs"] isa Vector{Observables.POVM}
    end

    @testset "Domain Errors" begin
        scenario = BipartiteNonSignaling(2,2,2,2)
        scenarioA = BipartiteNonSignaling(3,2,3,2)
        scenarioB = BipartiteNonSignaling(2,3,2,3)
        game = BellGame([0 1 1 0;0 0 0 1;0 0 0 1;1 0 0 1],2)
        ρ_AB = States.bell_states[1]
        POVMs = [
            Observables.POVM([[1 0;0 0],[0 0;0 1]]),
            Observables.POVM([[0.5 0.5;0.5 0.5],[0.5 -0.5;-0.5 0.5]])
        ]

        @test_throws DomainError Nonlocality.optimize_measurement(scenario, game, ρ_AB)

        @test_throws DomainError Nonlocality.optimize_measurement(scenarioA, game, ρ_AB, A_POVMs=POVMs)
        @test_throws DomainError Nonlocality.optimize_measurement(scenarioB, game, ρ_AB, B_POVMs=POVMs)
    end
end

end
