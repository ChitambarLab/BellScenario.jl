using Test, XPORTA, Polyhedra

@testset "./src/ConvexPolytope/adjacency_decomposition.jl" begin

using BellComm: ConvexPolytope

@testset "ConvexPolytope.rotate_facet()" begin

    F = HalfSpace([-1,0,0], 0)
    G1 = HalfSpace([0,0,-1], 0)
    G2 = HalfSpace([0,0,1], 1)

    xbar = [1,0,0]


    # neighboring facet 1
    @test ConvexPolytope.rotate_facet(F,G1,xbar) == HalfSpace([0,0,-1], 0)
    @test ConvexPolytope.rotate_facet(F,G2,xbar) == HalfSpace([0,-1,0], 0)


end

end
