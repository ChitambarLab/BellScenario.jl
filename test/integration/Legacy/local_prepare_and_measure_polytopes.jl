using Test, LinearAlgebra

@testset "test/integration/local_polytope_facets.jl" begin

using BellScenario: LocalPolytope, ConvexPolytope, Degeneracy, Behavior

dir = "./test/integration/files"

@testset "21->12 polytope" begin
    α = (2,1)
    β = (1,2)

    vertices = LocalPolytope.vertices(α, β, bits=1)

    # There should be (1*2)(1 + 1*(2 - 1)) vertices for fixed direction
    @test length(vertices) == 4
    @test size(vertices[1]) == (2,1)

    constraints = ConvexPolytope.facet_constraints(vertices, dir=dir)

    @test length(constraints[1]) == 4
    @test constraints[1] == [[0 -1 0],[0 0 -1],[-1 0 1],[-1 1 0]]
    @test constraints[2] == []
    @test constraints[3] == []

    canonical_groups = Degeneracy.canonical_facets(α,β,constraints[1],"normalized")

    @test issetequal(constraints[1], collect(Iterators.flatten(canonical_groups)))
    @test length(canonical_groups) == 1
    @test issetequal(canonical_groups[1], [
        [0 -1 0],[0 0 -1],[-1 0 1],[-1 1 0]
    ])
end

@testset "31->12 polytope" begin
    α = (3,1)
    β = (1,2)

    vertices = LocalPolytope.vertices(α, β, bits=1)

    # There should be (1*2)(1 + 1*(2 - 1)) vertices for fixed direction
    @test length(vertices) == 8
    @test size(vertices[1]) == (3,1)

    constraints = ConvexPolytope.facet_constraints(vertices, dir=dir)

    @test length(constraints[1]) == 6
    @test constraints[2] == []
    @test constraints[3] == []

    canonical_groups = Degeneracy.canonical_facets(α,β,constraints[1],"normalized")

    @test issetequal(constraints[1], collect(Iterators.flatten(canonical_groups)))
    @test length(canonical_groups) == 1
    @test issetequal(canonical_groups[1], [
        [0 -1 0 0],[0 0 -1 0],[0 0 0 -1],[-1 0 1 0],[-1 1 0 0],[-1 0 0 1]
    ])
end

@testset "31->13 polytope" begin
    α = (3,1)
    β = (1,3)

    vertices = LocalPolytope.vertices(α, β, bits=1)

    # There should be (1*3)(1 + 3*(3 - 1)) vertices for fixed direction
    @test length(vertices) == 21
    @test size(vertices[1]) == (6,1)

    constraints = ConvexPolytope.facet_constraints(vertices, dir=dir)

    @test length(constraints[1]) == 15
    @test constraints[2] == []
    @test constraints[3] == []

    canonical_groups = Degeneracy.canonical_facets(α,β,constraints[1],"normalized")

    @test issetequal(constraints[1], collect(Iterators.flatten(canonical_groups)))
    @test length(canonical_groups) == 2
    @test issetequal(canonical_groups[1], [
        [0 -1 0 0 0 0 0], [0 0 -1 0 0 0 0], [0 0 0 -1 0 0 0],
        [0 0 0 0 -1 0 0], [0 0 0 0 0 -1 0], [0 0 0 0 0 0 -1],
        [-1 0 0 1 0 0 1], [-1 0 1 0 0 1 0], [-1 1 0 0 1 0 0]
    ])
    @test issetequal(canonical_groups[2], [
        [-1 -1 0 1 -1 1 0], [-1 -1 1 0 -1 0 1], [-1 0 -1 1 1 -1 0],
        [-1 0 1 -1 1 0 -1], [-1 1 -1 0 0 -1 1], [-1 1 0 -1 0 1 -1]
    ])
end

@testset "41->14 polytope" begin
    α = (4,1)
    β = (1,4)

    vertices = LocalPolytope.vertices(α, β, bits=1)

    # There should be (1*4)(1 + 7*(4 - 1)) vertices for fixed direction
    @test length(vertices) == 88
    @test size(vertices[1]) == (12,1)

    constraints = ConvexPolytope.facet_constraints(vertices, dir=dir)

    @test length(constraints[1]) == 664
    @test constraints[2] == []
    @test constraints[3] == []

    canonical_groups = Degeneracy.canonical_facets(α,β,constraints[1],"normalized")

    ggame_facet = [-2 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1]*Behavior.norm_to_gen_proj(α,β)
    @test in(ggame_facet, canonical_groups[2])

    @test issetequal(constraints[1], collect(Iterators.flatten(canonical_groups)))
    @test length(canonical_groups) == 6
end

@testset "22->12 polytope" begin
    α = (2,2)
    β = (1,2)

    vertices = LocalPolytope.vertices(α, β, bits=1)

    # There should be (4*2)(1 + 1*(2 - 1)) vertices for fixed direction
    @test length(vertices) == 16
    @test size(vertices[1]) == (6,1)

    constraints = ConvexPolytope.facet_constraints(vertices, dir=dir)

    @test length(constraints[1]) == 8
    @test constraints[2] == []
    @test constraints[3] == []

    canonical_groups = Degeneracy.canonical_facets(α,β,constraints[1],"fixed-direction")

    @test issetequal(constraints[1], collect(Iterators.flatten(canonical_groups)))
    @test length(canonical_groups) == 1
end

@testset "22->22 polytope (fixed-direction)" begin
    α = (2,2)
    β = (2,2)

    vertices = LocalPolytope.vertices(α, β, bits=1, rep="fixed-direction")

    # There should be (4*4)(1 + 1*(4 - 1)) vertices for fixed direction
    @test length(vertices) == 64
    @test size(vertices[1]) == (10,1)

    constraints = ConvexPolytope.facet_constraints(vertices, dir=dir)
    @test constraints[2] == []
    @test constraints[3] == []

    @test length(constraints[1]) == 16

    canonical_groups = Degeneracy.canonical_facets(α,β,constraints[1],"fixed-direction")

    @test issetequal(constraints[1], collect(Iterators.flatten(canonical_groups)))
    @test length(canonical_groups) == 1
end

@testset "33->13 polytope" begin
    α = (3,3)
    β = (1,3)

    vertices = LocalPolytope.vertices(α, β, bits=1)

    # There should be (27*3)(1 + 3*(3 - 1)) vertices for fixed direction
    @test length(vertices) == 567
    @test size(vertices[1]) == (24,1)

    constraints = ConvexPolytope.facet_constraints(vertices, dir=dir)

    @test length(constraints[1]) == 33
    @test constraints[2] == []
    @test constraints[3] == []

    canonical_groups = Degeneracy.canonical_facets(α,β,constraints[1],"fixed-direction")

    @test issetequal(constraints[1], collect(Iterators.flatten(canonical_groups)))
    @test length(canonical_groups) == 2
end

@testset "32->22 polytope" begin
    α = (3,2)
    β = (2,2)

    vertices = LocalPolytope.vertices(α, β, bits=1)

    # There should be (27*3)(1 + 3*(3 - 1)) vertices for fixed direction
    @test length(vertices) == 320
    @test size(vertices[1]) == (15,1)

    constraints = ConvexPolytope.facet_constraints(vertices, dir=dir)

    @test length(constraints[1]) == 864
    @test constraints[2] == []
    @test constraints[3] == []

    canonical_groups = Degeneracy.canonical_facets(α,β,constraints[1],"fixed-direction")

    @test issetequal(constraints[1], collect(Iterators.flatten(canonical_groups)))
    @test length(canonical_groups) == 9
end

@testset "32->12 polytope" begin
    α = (3,2)
    β = (1,2)

    vertices = LocalPolytope.vertices(α, β, bits=1)

    @test length(vertices) == 64
    @test size(vertices[1]) == (9,1)

    constraints = ConvexPolytope.facet_constraints(vertices, dir=dir)

    @test length(constraints[1]) == 12
    @test constraints[2] == []
    @test constraints[3] == []

    canonical_groups = Degeneracy.canonical_facets(α,β,constraints[1],"fixed-direction")

    @test issetequal(constraints[1], collect(Iterators.flatten(canonical_groups)))
    @test length(canonical_groups) == 1
end

@testset "41->12 polytope" begin
    α = (4,1)
    β = (1,2)

    vertices = LocalPolytope.vertices(α, β, bits=1)

    @test length(vertices) == 16
    @test size(vertices[1]) == (4,1)

    constraints = ConvexPolytope.facet_constraints(vertices, dir=dir)

    @test length(constraints[1]) == 8
    @test constraints[2] == []
    @test constraints[3] == []

    canonical_groups = Degeneracy.canonical_facets(α,β,constraints[1],"normalized")

    @test issetequal(constraints[1], collect(Iterators.flatten(canonical_groups)))
    @test length(canonical_groups) == 1
end

@testset "33->12 polytope" begin
    α = (3,3)
    β = (1,2)

    vertices = LocalPolytope.vertices(α, β, bits=1)

    @test length(vertices) == 216
    @test size(vertices[1]) == (15,1)

    constraints = ConvexPolytope.facet_constraints(vertices, dir=dir)

    @test length(constraints[1]) == 18
    @test constraints[2] == []
    @test constraints[3] == []

    canonical_groups = Degeneracy.canonical_facets(α,β,constraints[1],"fixed-direction")

    @test issetequal(constraints[1], collect(Iterators.flatten(canonical_groups)))
    @test length(canonical_groups) == 1
end

@testset "32->13 polytope" begin
    α = (3,2)
    β = (1,3)

    vertices = LocalPolytope.vertices(α, β, bits=1)

    @test length(vertices) == 168
    @test size(vertices[1]) == (15,1)

    constraints = ConvexPolytope.facet_constraints(vertices, dir=dir)

    @test length(constraints[1]) == 24
    @test constraints[2] == []
    @test constraints[3] == []

    canonical_groups = Degeneracy.canonical_facets(α,β,constraints[1],"fixed-direction")

    @test issetequal(constraints[1], collect(Iterators.flatten(canonical_groups)))
    @test length(canonical_groups) == 2
end

@testset "41->13 polytope" begin
    α = (4,1)
    β = (1,3)

    vertices = LocalPolytope.vertices(α, β, bits=1)

    @test length(vertices) == 45
    @test size(vertices[1]) == (8,1)

    constraints = ConvexPolytope.facet_constraints(vertices, dir=dir)

    @test length(constraints[1]) == 36
    @test constraints[2] == []
    @test constraints[3] == []

    canonical_groups = Degeneracy.canonical_facets(α,β,constraints[1],"normalized")

    @test issetequal(constraints[1], collect(Iterators.flatten(canonical_groups)))
    @test length(canonical_groups) == 2
end

@testset "31->14 polytope" begin
    α = (3,1)
    β = (1,4)

    vertices = LocalPolytope.vertices(α, β, bits=1)

    @test length(vertices) == 40
    @test size(vertices[1]) == (9,1)

    constraints = ConvexPolytope.facet_constraints(vertices, dir=dir)

    @test length(constraints[1]) == 72
    @test constraints[2] == []
    @test constraints[3] == []

    canonical_groups = Degeneracy.canonical_facets(α,β,constraints[1],"normalized")

    @test issetequal(constraints[1], collect(Iterators.flatten(canonical_groups)))
    @test length(canonical_groups) == 3
end

end
