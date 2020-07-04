using Test, XPORTA

using Polyhedra: HalfSpace, vrep, convexhull, points

@testset "./src/ConvexPolytope/adjacency_decomposition.jl" begin

using BellScenario: ConvexPolytope

@testset "ConvexPolytope.rotate_facet()" begin
    @testset "simplex rotation" begin
        F = HalfSpace([-1,0,0], 0)
        G1 = HalfSpace([0,0,-1], 0)
        G2 = HalfSpace([0,0,1], 1)

        xbar = [1,0,0]

        # neighboring facet 1
        @test ConvexPolytope.rotate_facet(F,G1,xbar) == HalfSpace([0,0,-1], 0)
        @test ConvexPolytope.rotate_facet(F,G2,xbar) == HalfSpace([1,0,1], 1)
    end


    @testset "CHSH polytope" begin
        chsh_vertices = [
            -1 -1 -1 -1  1  1  1  1;
            -1 -1 -1  1  1 -1  1 -1;
            -1 -1  1 -1 -1  1 -1  1;
            -1 -1  1  1 -1 -1 -1 -1;
            -1  1 -1 -1  1  1 -1 -1;
            -1  1 -1  1  1 -1 -1  1;
            -1  1  1 -1 -1  1  1 -1;
            -1  1  1  1 -1 -1  1  1;
             1 -1 -1 -1 -1 -1  1  1;
             1 -1 -1  1 -1  1  1 -1;
             1 -1  1 -1  1 -1 -1  1;
             1 -1  1  1  1  1 -1 -1;
             1  1 -1 -1 -1 -1 -1 -1;
             1  1 -1  1 -1  1 -1  1;
             1  1  1 -1  1 -1  1 -1;
             1  1  1  1  1  1  1  1;
        ]

        chsh_inequalities = [
            -1 0 0 1 0 1 0 0 1;
            -1 0 1 0 1 0 0 0 1;
            0 -1 0 1 0 0 0 1 1;
            0 -1 1 0 0 0 1 0 1;
            0 1 -1 0 0 0 1 0 1;
            0 1 0 -1 0 0 0 1 1;
            0 1 0 1 0 0 0 -1 1;
            0 1 1 0 0 0 -1 0 1;
            1 0 -1 0 1 0 0 0 1;
            1 0 0 -1 0 1 0 0 1;
            1 0 0 1 0 -1 0 0 1;
            1 0 1 0 -1 0 0 0 1;
            -1 0 -1 0 -1 0 0 0 1;
            -1 0 0 -1 0 -1 0 0 1;
            0 -1 -1 0 0 0 -1 0 1;
            0 -1 0 -1 0 0 0 -1 1;
            0 0 0 0 -1 1 1 1 2;
            0 0 0 0 1 -1 1 1 2;
            0 0 0 0 1 1 -1 1 2;
            0 0 0 0 1 1 1 -1 2;
            0 0 0 0 -1 -1 -1 1 2;
            0 0 0 0 -1 -1 1 -1 2;
            0 0 0 0 -1 1 -1 -1 2;
            0 0 0 0 1 -1 -1 -1 2
        ]


        ch_vertices = BellScenario.LocalPolytope.vertices((2,2),(2,2))
        ch_ieq = traf(POI(vertices = vcat(map(v -> v', ch_vertices)...)))
        ch_ineqs = convert.(Int,ch_ieq.inequalities)

        ch_facets = [ HalfSpace(ch_ineqs[i,1:(end-1)], ch_ineqs[i,end]) for i in 1:size(ch_ineqs,1) ]



        F = HalfSpace(chsh_inequalities[1,1:(end-1)], chsh_inequalities[1,end])

        map( v -> F.a'*v, eachrow(chsh_vertices))

        verts = [chsh_vertices[i,:] for i in 1:size(chsh_vertices,1)]
        vert_rep = vrep(verts)

        filter(v -> F.a' * v == F.β, verts)

        convexhull(collect(eachrow(chsh_vertices))...)

        println(chsh_inequalities)

        ConvexPolytope.adjacent_facets(F, verts)

        hcat(verts...)'
    end

end

end
