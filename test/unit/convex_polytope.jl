using Test

@testset "/src/ConvexPolytope.jl" begin

using BellScenario: ConvexPolytope

@testset "linear_interpolation()" begin

    @test isapprox(ConvexPolytope.linear_interpolation([1;0;0],[0;1;0],dw=0.2), map(
        x -> reshape(x,(3,1)),
        [[1.0;0;0],[0.8;0.2;0],[0.6;0.4;0],[0.4;0.6;0],[0.2;0.8;0],[0;1.0;0]]
    ))

    @test all(map(
        x -> ConvexPolytope.facet_contains([0 0 0 -1], x),
        ConvexPolytope.linear_interpolation([1;0;0],[0;1;0],dw=0.2)
    ))
end

@testset "facet_contains()" begin

    @test ConvexPolytope.facet_contains([-1 1 1 1],[1/3;1/3;1/3])
    @test !(ConvexPolytope.facet_contains([-1 1 1 1],[1;1;1]))

    @test ConvexPolytope.facet_contains([-1 1 0 -1 0 1 -1],[2/3;1/6;1/6;1/6;2/3;1/6])
end

@testset "violated_facets()" begin

    facets = [[0 -1 0 0],[0 0 -1 0],[0 0 0 -1],[-1 1 1 1]]
    @test ConvexPolytope.violated_facets(facets, [1;1;1]) == [[-1 1 1 1]]
    @test ConvexPolytope.violated_facets(facets, [1/4;1/4;1/4]) == []

    facets = [
        [0 -1 0 0 0 0 0],[0 0 -1 0 0 0 0],[0 0 0 -1 0 0 0],[0 0 0 0 -1 0 0],
        [0 0 0 0 0 -1 0],[0 0 0 0 0 0 -1],[-1 0 0 1 0 0 1],[-1 0 1 0 0 1 0],
        [-1 1 0 0 1 0 0],[-1 -1 0 1 -1 1 0],[-1 -1 1 0 -1 0 1],[-1 0 -1 1 1 -1 0],
        [-1 0 1 -1 1 0 -1], [-1 1 -1 0 0 -1 1], [-1 1 0 -1 0 1 -1]
    ]
    @test ConvexPolytope.violated_facets(facets, [1;0;0;0;1;0]) == [[-1 1 0 -1 0 1 -1]]
end

@testset "facet_violation()" begin

    @test ConvexPolytope.facet_violation([-1 1 0 -1 0 1 -1], [1;1;0;0;0;1;0]) == 1
    @test ConvexPolytope.facet_violation([-1 1 0 -1 0 1 -1], [1;0;0;0;1;0]) == 1

    @test_throws ArgumentError ConvexPolytope.facet_violation([-1 1 0 -1 0 1 -1],[1;1])

    @test ConvexPolytope.facet_violation([-1 1 1 1], [1;1;1]) == 2
end

@testset "facet_distance()" begin
    @test ConvexPolytope.facet_distance([0 -1 0 0],[1;1;1]) ≈ -1
    @test ConvexPolytope.facet_distance([-1 1 1 1], [1;1;1]) ≈ 2/sqrt(3)
end

@testset "contains()" begin

    vertices = [[0;0;0],[1;0;0],[0;1;0],[0;0;1]]
    facets = [[0 -1 0 0],[0 0 -1 0],[0 0 0 -1],[-1 1 1 1]]

    @test !(ConvexPolytope.facets_contain(facets, [1;1;1]))
    @test ConvexPolytope.facets_contain(facets,[1/4;1/4;1/4])
    @test ConvexPolytope.facets_contain(facets,[1/3;1/3;1/3])
    @test all(map(
        (v) -> ConvexPolytope.facets_contain(facets, v),
        vertices
    ))

end

@testset "collect_vertices()" begin
    facets = [[-1 1 0 -1 0 1 -1],[0 -1 0 0 0 0 0]]
    vertices = [
        [1;0;0;0;0;0;0],[1;1;1;1;0;0;0],[1;0;0;0;1;1;1],[1;1;1;0;0;0;0],
        [1;1;1;0;0;0;1],[1;1;0;1;0;0;0],[1;1;0;1;0;1;0],[1;0;1;1;0;0;0],
        [1;0;1;1;1;0;0]
    ]

    facet_dict = ConvexPolytope.collect_vertices(facets, vertices)

    @test length(facet_dict) == 2
    @test facet_dict[1]["facet"] == [-1 1 0 -1 0 1 -1]
    @test facet_dict[1]["vertices"] == [[1;1;1;0;0;0;0],[1;1;0;1;0;1;0]]
    @test facet_dict[2]["facet"] == [0 -1 0 0 0 0 0]
    @test facet_dict[2]["vertices"] == [[1;0;0;0;0;0;0],[1;0;0;0;1;1;1],[1;0;1;1;0;0;0],[1;0;1;1;1;0;0]]
end

@testset "collect_facets()" begin
    facets = [[-1 1 0 -1 0 1 -1],[0 -1 0 0 0 0 0]]
    vertices = [
        [1;0;0;0;0;0;0],[1;0;0;0;1;1;1],[1;1;1;0;0;0;1],[1;0;0;0;1;1;0]
    ]

    facet_dict = ConvexPolytope.collect_facets(vertices, facets)

    @test length(facet_dict) == 4
    @test facet_dict[1]["vertex"] == [1;0;0;0;0;0;0]
    @test facet_dict[1]["facets"] == [[0 -1 0 0 0 0 0]]
    @test facet_dict[2]["vertex"] == [1;0;0;0;1;1;1]
    @test facet_dict[2]["facets"] == [[0 -1 0 0 0 0 0]]
    @test facet_dict[3]["vertex"] == [1;1;1;0;0;0;1]
    @test facet_dict[3]["facets"] == []
    @test facet_dict[4]["vertex"] == [1;0;0;0;1;1;0]
    @test facet_dict[4]["facets"] == [[-1 1 0 -1 0 1 -1],[0 -1 0 0 0 0 0]]
end

@testset "contains()" begin
    facets = [[0 -1 -1 -1 0 0 0], [0 0 0 0 -1 -1 -1], [0 0 0 0 0 0 0]]

    @test ConvexPolytope.contains(facets, [1;0;0;0;0;0;0])
    @test ConvexPolytope.contains(facets, [1;1;1;1;0;0;0])
    @test ConvexPolytope.contains(facets, [1;0.5;0.5;0.5;0.5;0.5;0.5])
    @test !ConvexPolytope.contains(facets, [1;-1;0;0;0;0;0])
end

end
