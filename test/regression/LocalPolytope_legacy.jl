using Test, LinearAlgebra

@testset "./src/LocalPolytope/vertices.jl - Legacy Regression" begin

using BellScenario

@testset "vertices(scenario :: BipartiteNoSignaling)" begin

    @testset "no-signaling rep" begin
        nosignaling_reshape_old_vertices(scenario, v_old) = begin
            dim_α = (scenario.A-1)*(scenario.X)
            dim_β = (scenario.B-1)*(scenario.Y)
            dim_αβ = dim_α*dim_β
            dim_ns = dim_α + dim_β + dim_αβ

            map( v -> begin
                v_reshape = zeros(Int64, dim_ns)
                v_reshape[1:dim_α] = reshape(v[1:dim_α],(scenario.X,scenario.A-1))'[:]
                v_reshape[1+dim_α:dim_α+dim_β] = reshape(v[1+dim_α:dim_α+dim_β],(scenario.Y,scenario.B-1))'[:]
                v_reshape[1+dim_α+dim_β:dim_ns] = reshape(v[1+dim_α+dim_β:dim_ns], (scenario.X*scenario.Y,(scenario.B-1)*(scenario.A-1)))'[:]

                v_reshape
            end, v_old)
        end

        @testset "CHSH case" begin
            scenario = BipartiteNoSignaling(2,2,2,2)
            v_new = LocalPolytope.vertices(scenario, "no-signaling")
            v_old = nosignaling_reshape_old_vertices(scenario, LocalPolytope.vertices((2,2),(2,2),dits=1,rep="no-signaling"))

            @test sort(v_new) == sort(v_old)
        end

        @testset "(4,3)-(4,3) case" begin
            scenario = BipartiteNoSignaling(3,3,4,4)
            v_new = LocalPolytope.vertices(scenario, "no-signaling")

            v_old = nosignaling_reshape_old_vertices(scenario, LocalPolytope.vertices((4,3),(4,3),dits=1,rep="no-signaling"))

            @test sort(v_new) == sort(v_old)
        end

        @testset "(4,4)-(3,3) case" begin
            scenario = BipartiteNoSignaling(4,3,4,3)
            v_new = LocalPolytope.vertices(scenario, "no-signaling")

            v_old = nosignaling_reshape_old_vertices(scenario, LocalPolytope.vertices((4,4),(3,3),dits=1,rep="no-signaling"))

            @test sort(v_new) == sort(v_old)
        end

    end

    @testset "normalized rep" begin
        normalized_reshape_old_vertices(scenario, v_old) = begin
            map(v -> reshape(convert.(Int64, v), (scenario.X*scenario.Y, scenario.A*scenario.B-1))'[:], v_old)
        end

        @testset "CHSH" begin
            scenario = BipartiteNoSignaling(2,2,2,2)
            v_new = LocalPolytope.vertices(scenario, "normalized")
            v_old = normalized_reshape_old_vertices(scenario, LocalPolytope.vertices((2,2),(2,2),dits=1,rep="normalized"))

            @test sort(v_new) == sort(v_old)
        end

        @testset "(4,3)-(4,3) case" begin
            scenario = BipartiteNoSignaling(3,3,4,4)
            v_new = LocalPolytope.vertices(scenario, "normalized")
            v_old = normalized_reshape_old_vertices(scenario, LocalPolytope.vertices((4,3),(4,3),dits=1,rep="normalized"))

            @test sort(v_new) == sort(v_old)
        end

        @testset "(4,4)-(3,3) case" begin
            scenario = BipartiteNoSignaling(4,3,4,3)
            v_new = LocalPolytope.vertices(scenario, "normalized")
            v_old = normalized_reshape_old_vertices(scenario, LocalPolytope.vertices((4,4),(3,3),dits=1,rep="normalized"))

            @test sort(v_new) == sort(v_old)
        end
    end

    @testset "generalized rep" begin
        generalized_reshape_old_vertices(scenario, v_old) = begin
            map(v -> reshape(convert.(Int64, v), (scenario.X*scenario.Y, scenario.A*scenario.B))'[:], v_old)
        end

        @testset "CHSH case" begin
            scenario = BipartiteNoSignaling(2,2,2,2)
            v_new = LocalPolytope.vertices(scenario, "generalized")
            v_old = generalized_reshape_old_vertices(scenario, LocalPolytope.vertices((2,2),(2,2),dits=1,rep="generalized"))

            @test sort(v_new) == sort(v_old)
        end

        @testset "(4,3)-(4,3) case" begin
            scenario = BipartiteNoSignaling(3,3,4,4)
            v_new = LocalPolytope.vertices(scenario, "generalized")
            v_old = generalized_reshape_old_vertices(scenario, LocalPolytope.vertices((4,3),(4,3),dits=1,rep="generalized"))

            @test sort(v_new) == sort(v_old)
        end

        @testset "(4,4)-(3,3) case" begin
            scenario = BipartiteNoSignaling(4,3,4,3)
            v_new = LocalPolytope.vertices(scenario, "generalized")
            v_old = generalized_reshape_old_vertices(scenario, LocalPolytope.vertices((4,4),(3,3),dits=1,rep="generalized"))

            @test sort(v_new) == sort(v_old)
        end
    end
end

end
