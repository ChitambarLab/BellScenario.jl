using Test, QBase

@testset "./src/Nonlocality/optimize_measurements.jl" begin

using BellScenario

@testset "Nonlocality.optimize_measurement(LocalSignaling)" begin
    @testset "non State types" begin
        scenario = LocalSignaling(3,3,2)
        game = BellGame([1 0 0;0 1 0;0 0 1], 2)
        ρ_states = [[1 0;0 0],[0 0;0 1],[0 0;0 1]]

        dict = Nonlocality.optimize_measurement(scenario, game, ρ_states)

        @test isapprox(dict["violation"], 0.0, atol=1e-4)
        @test isapprox(dict["povm"], [[1 0;0 0],[0 0;0 0.5],[0 0;0 0.5]], atol=1e-4)
    end

    @testset "trine states" begin
        scenario = LocalSignaling(3,3,2)
        game = BellGame([1 0 0;0 1 0;0 0 1], 2)
        ρ_states = trine_qubit_states()

        dict = Nonlocality.optimize_measurement(scenario, game, ρ_states)

        @test isapprox(dict["violation"], 0.0, atol=1e-5)
        @test all(isapprox.(dict["povm"][1], 2/3*ρ_states[1], atol=1e-5))
        @test all(isapprox.(dict["povm"][2], 2/3*ρ_states[2], atol=1e-5))
        @test all(isapprox.(dict["povm"][3], 2/3*ρ_states[3], atol=1e-5))

        @test dict["states"] == ρ_states
        @test dict["game"] == game
        @test dict["scenario"] == scenario
    end

    @testset "bb84 states" begin
        scenario = LocalSignaling(4,4,2)
        game = BellGame([1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1], 2)
        ρ_states = bb84_qubit_states()

        dict = Nonlocality.optimize_measurement(scenario, game, ρ_states)

        @test isapprox(dict["violation"], 0.0, atol=1e-5)
        @test all(isapprox.(dict["povm"][1], 1/2*ρ_states[1], atol=1e-3))
        @test all(isapprox.(dict["povm"][2], 1/2*ρ_states[2], atol=1e-3))
        @test all(isapprox.(dict["povm"][3], 1/2*ρ_states[3], atol=1e-3))
        @test all(isapprox.(dict["povm"][4], 1/2*ρ_states[4], atol=1e-3))
    end

    @testset "sic qubit states" begin
        scenario = LocalSignaling(4,4,2)
        game = BellGame([1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1], 2)
        ρ_states = sic_qubit_states()

        dict = Nonlocality.optimize_measurement(scenario, game, ρ_states)

        @test isapprox(dict["violation"], 0.0, atol=1e-5)
        @test all(isapprox.(dict["povm"][1], 1/2*ρ_states[1], atol=1e-5))
        @test all(isapprox.(dict["povm"][2], 1/2*ρ_states[2], atol=1e-5))
        @test all(isapprox.(dict["povm"][3], 1/2*ρ_states[3], atol=1e-5))
        @test all(isapprox.(dict["povm"][4], 1/2*ρ_states[4], atol=1e-5))
    end

    @testset "Errors" begin
        scenario = LocalSignaling(3,3,2)
        states = bb84_qubit_states()
        game = BellGame([1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1], 2)

        @test_throws DomainError Nonlocality.optimize_measurement(scenario, game, states)

        scenario = LocalSignaling(4,4,3)
        @test_throws DomainError  Nonlocality.optimize_measurement(scenario, game, states)

        scenario = LocalSignaling(3,3,2)
        game = BellGame([1 0 0;0 1 0;0 0 1],2)
        states = [[1 0;0 1],[0 0;0 1],[0.5 0;0 0.5]]
        @test_throws DomainError Nonlocality.optimize_measurement(scenario, game, states)
    end
end

@testset "Nonlocality.optimize_measurement(::BipartiteNonSignaling)" begin
    @testset "matrix states and povms" begin
        scenario = BipartiteNonSignaling(2,2,2,2)
        game = BellGame([0 1 1 0;0 0 0 1;0 0 0 1;1 0 0 1],2)
        ρ_AB = [1 0 0 1;0 0 0 0;0 0 0 0;1 0 0 1]/2
        POVMs = [
            [[1 0;0 0],[0 0;0 1]],
            [[0.5 0.5;0.5 0.5],[0.5 -0.5;-0.5 0.5]]
        ]

        opt_dictA = Nonlocality.optimize_measurement(scenario, game, ρ_AB, A_POVMs=POVMs)
        opt_dictB = Nonlocality.optimize_measurement(scenario, game, ρ_AB, B_POVMs=POVMs)

        @test isapprox(opt_dictA["score"], 2.20710279, atol=1e-5)
        @test isapprox(opt_dictA["violation"], 0.207102796, atol=1e-5)
        @test opt_dictA["A_POVMs"] isa Vector{<:POVM}
        @test opt_dictA["state"] isa State

        @test isapprox(opt_dictB["score"], 2.20710279, atol=1e-5)
        @test isapprox(opt_dictB["violation"], 0.207102796, atol=1e-5)
        @test opt_dictB["B_POVMs"] isa Vector{<:POVM}
        @test opt_dictB["state"] isa State
    end

    @testset "CHSH inequality" begin
        scenario = BipartiteNonSignaling(2,2,2,2)
        game = BellGame([0 1 1 0;0 0 0 1;0 0 0 1;1 0 0 1],2)
        ρ_AB = bell_states()[1]
        POVMs = [
            POVM([[1 0;0 0],[0 0;0 1]]),
            POVM([[0.5 0.5;0.5 0.5],[0.5 -0.5;-0.5 0.5]])
        ]

        @test game.β == 2

        opt_dictA = Nonlocality.optimize_measurement(scenario, game, ρ_AB, A_POVMs=POVMs)
        opt_dictB = Nonlocality.optimize_measurement(scenario, game, ρ_AB, B_POVMs=POVMs)

        @test length(keys(opt_dictA)) == 7
        @test opt_dictA["scenario"] == scenario
        @test opt_dictA["game"] == game
        @test opt_dictA["game"].β == 2
        @test isapprox(opt_dictA["score"], 2.20710279, atol=1e-5)
        @test isapprox(opt_dictA["violation"], 0.207102796, atol=1e-5)
        @test opt_dictA["state"] == ρ_AB
        @test opt_dictA["A_POVMs"] == POVMs
        @test isapprox(
            opt_dictA["B_POVMs"],
            [
                [bloch_qubit_state(π/4,0), bloch_qubit_state(5π/4,0)],
                [bloch_qubit_state(7π/4,0), bloch_qubit_state(3π/4,0)]
            ],
            atol=1e-5
        )

        @test length(keys(opt_dictB)) == 7
        @test opt_dictB["scenario"] == scenario
        @test opt_dictB["game"] == game
        @test opt_dictB["game"].β == 2
        @test isapprox(opt_dictB["score"], 2.2071028, atol=1e-5)
        @test isapprox(opt_dictB["violation"], 0.207102796, atol=1e-5)
        @test opt_dictB["state"] == ρ_AB
        @test opt_dictB["B_POVMs"] == POVMs
        @test isapprox(
            opt_dictB["A_POVMs"],
            [
                [bloch_qubit_state(π/4,0), bloch_qubit_state(5π/4,0)],
                [bloch_qubit_state(7π/4,0),bloch_qubit_state(3π/4,0)]
            ],
            atol=1e-5
        )
    end

    @testset "Domain Errors" begin
        scenario = BipartiteNonSignaling(2,2,2,2)
        scenarioA = BipartiteNonSignaling(3,2,3,2)
        scenarioB = BipartiteNonSignaling(2,3,2,3)
        game = BellGame([0 1 1 0;0 0 0 1;0 0 0 1;1 0 0 1],2)
        ρ_AB = bell_states()[1]
        POVMs = [
            POVM([[1 0;0 0],[0 0;0 1]]),
            POVM([[0.5 0.5;0.5 0.5],[0.5 -0.5;-0.5 0.5]])
        ]

        @test_throws DomainError Nonlocality.optimize_measurement(scenario, game, ρ_AB)

        @test_throws DomainError Nonlocality.optimize_measurement(scenarioA, game, ρ_AB, A_POVMs=POVMs)
        @test_throws DomainError Nonlocality.optimize_measurement(scenarioB, game, ρ_AB, B_POVMs=POVMs)

        @test_throws DomainError Nonlocality.optimize_measurement(
            scenario, game, [1 0 0 1;0 0 0 0;0 0 0 0;1 0 0 1],A_POVMs=POVMs
        )
        @test_throws DomainError Nonlocality.optimize_measurement(
            scenario, game, ρ_AB, A_POVMs = [[[1 0;0 1],[1 0;0 1]],[[1 0;0 0],[1 0;0 0]]]
        )
        @test_throws DomainError Nonlocality.optimize_measurement(
            scenario, game, ρ_AB, B_POVMs = [[[1 0;0 1],[1 0;0 1]],[[1 0;0 0],[1 0;0 0]]]
        )
    end
end

end
