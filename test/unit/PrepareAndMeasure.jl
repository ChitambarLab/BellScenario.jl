using Test

@testset "/src/PrepareAndMeasure.jl" begin

using BellScenario: PrepareAndMeasure

@testset "num_zero_rows()" begin

    @test 0 == PrepareAndMeasure.num_zero_rows([1 0;0 1])
    @test 1 == PrepareAndMeasure.num_zero_rows([1 1;0 0])
    @test 2 == PrepareAndMeasure.num_zero_rows([0 0;0 0])

    @test 0 == PrepareAndMeasure.num_zero_rows([1 0 1;0 1 0])
    @test 1 == PrepareAndMeasure.num_zero_rows([0 0 0;1 1 1])

    @test 1 == PrepareAndMeasure.num_zero_rows([0 1;0 0;1 0])
    @test 3 == PrepareAndMeasure.num_zero_rows([1 0;0 0;0 1;0 0;0 0])

    @test 2 == PrepareAndMeasure.num_zero_rows([1 0 0]')
end

@testset "vertices_α_interpolable()" begin
    @test PrepareAndMeasure.vertices_α_interpolable([1 0 0;0 1 1;0 0 0], [1 0 0;0 1 1;0 0 0])
    @test PrepareAndMeasure.vertices_α_interpolable([1 0 0;0 1 1;0 0 0], [0 1 0;1 0 1;0 0 0])
    @test PrepareAndMeasure.vertices_α_interpolable([1 1 0;0 0 1;0 0 0], [0 0 0;1 1 1;0 0 0])
    @test PrepareAndMeasure.vertices_α_interpolable([0 0 0;1 1 1;0 0 0], [1 1 0;0 0 1;0 0 0])

    @test !PrepareAndMeasure.vertices_α_interpolable([1 0 0;0 1 1;0 0 0], [0 0 0;1 0 0;0 1 1])
    @test !PrepareAndMeasure.vertices_α_interpolable([1 0 1;0 1 0;0 0 0], [0 0 0;0 1 0;1 0 1])
end

@testset "vertices_β_interpolable()" begin
    @test PrepareAndMeasure.vertices_β_interpolable([1 0 0;0 1 1;0 0 0], [1 0 0;0 1 1;0 0 0])
    @test PrepareAndMeasure.vertices_β_interpolable([1 1 0;0 0 1;0 0 0], [0 0 0;1 1 1;0 0 0])
    @test PrepareAndMeasure.vertices_β_interpolable([0 0 0;1 1 1;0 0 0], [1 1 0;0 0 1;0 0 0])
    @test PrepareAndMeasure.vertices_β_interpolable([1 0 0;0 1 1;0 0 0], [0 0 0;1 0 0;0 1 1])
    @test PrepareAndMeasure.vertices_β_interpolable([1 0 1;0 1 0;0 0 0], [0 0 0;0 1 0;1 0 1])
    @test PrepareAndMeasure.vertices_β_interpolable([1 0 0;0 1 1;0 0 0], [0 0 0;0 0 0;1 1 1])

    @test !PrepareAndMeasure.vertices_β_interpolable([1 0 0;0 1 1;0 0 0], [0 1 0;1 0 1;0 0 0])
    @test !PrepareAndMeasure.vertices_β_interpolable([1 0 0;0 1 1;0 0 0], [0 1 0;0 0 0;1 0 1])
end

@testset "vertices_interpolable()" begin
    @test PrepareAndMeasure.vertices_interpolable([1 0 0;0 1 1;0 0 0], [0 1 0;1 0 1;0 0 0])
    @test PrepareAndMeasure.vertices_interpolable([1 0 0;0 1 1;0 0 0], [1 0 0;0 1 1;0 0 0])
    @test PrepareAndMeasure.vertices_interpolable([1 1 0;0 0 1;0 0 0], [0 0 0;1 1 1;0 0 0])

    @test !PrepareAndMeasure.vertices_interpolable([1 0 1;0 1 0;0 0 0], [1 0 0;0 0 0;0 1 1])
end

end
