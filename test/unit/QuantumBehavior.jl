using Test

@testset "./src/Scenario/QuantumBehavior.jl" begin

using BellScenario: QuantumBehavior

@testset "simple separable cases" begin

    test_cases = [
        ([1;0;0;0],[1;1;1;1;1;1;1;1;1]),
        ([0;1;0;0],[1;1;1;0;0;0;0;0;0]),
        ([0;0;1;0],[1;0;0;1;1;0;0;0;0]),
        ([0;0;0;1],[1;0;0;0;0;0;0;0;0])
    ]

    @testset "ψ = $(case[1])" for case in test_cases
        povm = [[1 0;0 0],[0 0;0 1]]

        behavior = QuantumBehavior.bipartite_scenario([povm,povm], [povm,povm], case[1])

        @test behavior == reshape(case[2],(9,1))
    end
end

@testset "prepare_and_measure()" begin
    @test [1 1 0]' == QuantumBehavior.prepare_and_measure(
        [[1 0;0 0],[0 0;0 1]],
        [[[1 0;0 0],[0 0;0 1]]]
    )
    @test [1 1 0.5 0 0.5]' == QuantumBehavior.prepare_and_measure(
        [[1 0;0 0],[0 0;0 1]],
        [[[1 0;0 0],[0 0;0 1]], [[0.5 0.5;0.5 0.5], [0.5 -0.5;-0.5 0.5]]]
    )
end

end
