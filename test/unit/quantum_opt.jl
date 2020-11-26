using Test, QBase

@testset "./src/quantum_opt.jl" begin

using BellScenario

@testset "optimize_measurement(LocalSignaling)" begin
    @testset "trine states" begin
        scenario = LocalSignaling(3,3,2)
        game = BellGame([1 0 0;0 1 0;0 0 1], 2)
        ρ_states = States.trine_qubits

        dict = optimize_measurement(game, ρ_states, scenario)

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

        dict = optimize_measurement(game, ρ_states, scenario)

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

        dict = optimize_measurement(game, ρ_states, scenario)

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

        @test_throws DomainError optimize_measurement(game, states, scenario)

        scenario = LocalSignaling(4,4,3)
        @test_throws DomainError  optimize_measurement(game, states, scenario)
    end
end

end
