using Test

@testset "./src/LocalPolytope/linear_nonclassicality_witnesses.jl" begin

using BellScenario

@testset "linear_nonclassicality_witness" begin
    @testset "CHSH Inequality Correlation Vertices" begin
        
        chsh_correlation_vertices = [
            [-1, -1, -1, -1,  1,  1,  1,  1],
            [-1, -1, -1,  1,  1, -1,  1, -1],
            [-1, -1,  1, -1, -1,  1, -1,  1],
            [-1, -1,  1,  1, -1, -1, -1, -1],
            [-1,  1, -1, -1,  1,  1, -1, -1],
            [-1,  1, -1,  1,  1, -1, -1,  1],
            [-1,  1,  1, -1, -1,  1,  1, -1],
            [-1,  1,  1,  1, -1, -1,  1,  1],
            [ 1, -1, -1, -1, -1, -1,  1,  1],
            [ 1, -1, -1,  1, -1,  1,  1, -1],
            [ 1, -1,  1, -1,  1, -1, -1,  1],
            [ 1, -1,  1,  1,  1,  1, -1, -1],
            [ 1,  1, -1, -1, -1, -1, -1, -1],
            [ 1,  1, -1,  1, -1,  1, -1,  1],
            [ 1,  1,  1, -1,  1, -1,  1, -1],
            [ 1,  1,  1,  1,  1,  1,  1,  1],
        ]

        pr_box_test_behavior = [0, 0, 0, 0, 1, 1, 1, -1]

        facet_vec = LocalPolytope.linear_nonclassicality_witness(chsh_correlation_vertices, pr_box_test_behavior[:])

        @test facet_vec ≈ [0, 0, 0, 0, 0.5, 0.5, 0.5, -0.5, 1]  # [A0, A1, B0, B1, AB00, AB01, AB10, AB11]
    end
    
    @testset "CH Inequality Probability Vertices" begin
        chsh_scenario = BipartiteNonSignaling(2, 2, 2, 2)

        chsh_vertices = LocalPolytope.vertices(chsh_scenario, "non-signaling")

        pr_box_test_behavior = [1, 1, 1, 1, 1, 1, 1, 0] / 2  # nonlocal test behavior
        local_test_behavior = [1/2, 1/2, 1/2, 1/2, 1/4, 1/4, 1/4, 1/4]  # white noise behavior

        chsh_facet_vec = LocalPolytope.linear_nonclassicality_witness(chsh_vertices, pr_box_test_behavior[:])
        local_game_vec = LocalPolytope.linear_nonclassicality_witness(chsh_vertices, local_test_behavior[:])

        # -2*PA(0|0) - 2*PB(0|0) + 2*PAB(00|00) + 2*PAB(01|00) + 2*PAB(10|00) - 2*PAB(11|00) <= 0
        @test chsh_facet_vec ≈ [-2, 0, -2, 0, 2, 2, 2, -2, 0]
        @test local_game_vec ≈ [0, 0, 0, 0, 0, 0, 0, 0, 0]
    end
end

end