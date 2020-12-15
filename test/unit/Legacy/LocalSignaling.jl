using Test

@testset "./src/Legacy/LocalSignaling.jl" begin

using BellScenario

@testset "num_zero_rows()" begin

    @test 0 == BellScenario.num_zero_rows([1 0;0 1])
    @test 1 == BellScenario.num_zero_rows([1 1;0 0])
    @test 2 == BellScenario.num_zero_rows([0 0;0 0])

    @test 0 == BellScenario.num_zero_rows([1 0 1;0 1 0])
    @test 1 == BellScenario.num_zero_rows([0 0 0;1 1 1])

    @test 1 == BellScenario.num_zero_rows([0 1;0 0;1 0])
    @test 3 == BellScenario.num_zero_rows([1 0;0 0;0 1;0 0;0 0])

    @test 2 == BellScenario.num_zero_rows([1 0 0]')
end

@testset "vertices_α_interpolable()" begin
    @test BellScenario.vertices_α_interpolable([1 0 0;0 1 1;0 0 0], [1 0 0;0 1 1;0 0 0])
    @test BellScenario.vertices_α_interpolable([1 0 0;0 1 1;0 0 0], [0 1 0;1 0 1;0 0 0])
    @test BellScenario.vertices_α_interpolable([1 1 0;0 0 1;0 0 0], [0 0 0;1 1 1;0 0 0])
    @test BellScenario.vertices_α_interpolable([0 0 0;1 1 1;0 0 0], [1 1 0;0 0 1;0 0 0])

    @test !BellScenario.vertices_α_interpolable([1 0 0;0 1 1;0 0 0], [0 0 0;1 0 0;0 1 1])
    @test !BellScenario.vertices_α_interpolable([1 0 1;0 1 0;0 0 0], [0 0 0;0 1 0;1 0 1])
end

@testset "vertices_β_interpolable()" begin
    @test BellScenario.vertices_β_interpolable([1 0 0;0 1 1;0 0 0], [1 0 0;0 1 1;0 0 0])
    @test BellScenario.vertices_β_interpolable([1 1 0;0 0 1;0 0 0], [0 0 0;1 1 1;0 0 0])
    @test BellScenario.vertices_β_interpolable([0 0 0;1 1 1;0 0 0], [1 1 0;0 0 1;0 0 0])
    @test BellScenario.vertices_β_interpolable([1 0 0;0 1 1;0 0 0], [0 0 0;1 0 0;0 1 1])
    @test BellScenario.vertices_β_interpolable([1 0 1;0 1 0;0 0 0], [0 0 0;0 1 0;1 0 1])
    @test BellScenario.vertices_β_interpolable([1 0 0;0 1 1;0 0 0], [0 0 0;0 0 0;1 1 1])

    @test !BellScenario.vertices_β_interpolable([1 0 0;0 1 1;0 0 0], [0 1 0;1 0 1;0 0 0])
    @test !BellScenario.vertices_β_interpolable([1 0 0;0 1 1;0 0 0], [0 1 0;0 0 0;1 0 1])
end

@testset "vertices_interpolable()" begin
    @test BellScenario.vertices_interpolable([1 0 0;0 1 1;0 0 0], [0 1 0;1 0 1;0 0 0])
    @test BellScenario.vertices_interpolable([1 0 0;0 1 1;0 0 0], [1 0 0;0 1 1;0 0 0])
    @test BellScenario.vertices_interpolable([1 1 0;0 0 1;0 0 0], [0 0 0;1 1 1;0 0 0])

    @test !BellScenario.vertices_interpolable([1 0 1;0 1 0;0 0 0], [1 0 0;0 0 0;0 1 1])
end

end
