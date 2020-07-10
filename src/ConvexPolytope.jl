"""
This BellComm submodule provides methods for polyhedral analysis and transformations.
"""
module ConvexPolytope

using XPORTA, LinearAlgebra

# """
# facet_constraints(vertices, filename)
#
#     Uses xporta software to compute the facet inequalities given a set of
#     polytope vertices.
#
# Side-effect:
#     filename.poi and filename.poi.ieq files are saved to porta_tmp. porta_tmp
#     is optionally deleted at the end of function's execution.
#
# Inputs:
#     vertices: Array, contains polytope vertices.
#     filename: String (optional), name of .poi and .ieq files
#     dir: String (optional), directory where porta_tmp is created
#     cleanup: Boolean, if true tmp files are removed
#     suppress: Boolean, if true xporta does not print to STDOUT
#
# Outputs:
#     constraints: Tuple, contains arrays of dual vectors to vertices, (e.g. facets).
#       The constraint tupe is ordered in the following way
#           1: less than constraints <= 0
#           2: equal to constraints == 0
#           3: greater than constraints >= 0
# """
function facet_constraints(vertices; filename="facet_constraints_tmp", dir="./", cleanup=true, suppress=true)
    ieq = traf(POI( vertices=convert.(Int, vcat(map(v -> v', vertices)...))), filename=filename, dir=dir, cleanup=cleanup, verbose=!suppress)

    equalities = map(eq -> cat(-1*eq[end], eq[1:end-1]', dims=2) , eachrow(ieq.equalities))
    inequalities = map(ineq -> cat(-1*ineq[end], ineq[1:end-1]', dims=2) , eachrow(ieq.inequalities))

    (inequalities, equalities, [])
end

# """
# linear_interpolation(x1, x2; dw=0.1)
#
#     Creates an array of points which lie on the line connecting x1 to x2. Points
#     are interpolated by a step interval of dw where 0:dw:1.
#
# Inputs:
#     x1: Col Vector, point in convex polytope space
#     x2: Col Vector, point in convex polytope space
#     dw: float, 0<dw<1, sets the number of steps used in the interpolation
#
# Output:
#     interpolation: Array, contains points in convex polytope which scan across
#                    the line starting from x1 going to x2.
# """
function linear_interpolation(x1, x2; dw=0.1)
    ws = 0:dw:1
    dim = length(x1)
    interpolation = map(w -> reshape((1 - w)*x1 + w*x2,(dim, 1)), ws)

    interpolation
end

# """
# facet_contains(facet, x)
#
#     Returns true if the point x lies approximately on the facet.
#
# Inputs:
#     facet: Row Vector, contains hyperplane constant and normal vector of hyperplane.
#     x: Col Vector, a point the convex polytope space.
#
# Outputs:
#     Boolean, true if point x lies on the facet.
# """
function facet_contains(facet, x)
    isapprox(1 + facet_violation(facet, x), 1)
end

# """
# contains(facets, x):
#
#     Checks if the convex polytope contains the point x by checking each facet
#     for violation.
#
# Inputs:
#     facets: Array, contains facet inequalities
#     x: Col Vector, the point to test
#
# Output:
#     Boolean, true if the polytope contains point x
# """
function contains(facets, x)
    length(violated_facets(facets, x)) == 0
end

# """
# violated_facets(facets, x)
#
#     Returns the list of facets which are violated by the provided point, x. If
#     no facets are violated, an empty array is returned.
#
# Inputs:
#     facet: Row Vector, contains hyperplane constant and normal vector of hyperplane.
#     x: Col Vector, a point the convex polytope space.
#
# Outputs:
#     violations: Array, contains the facet inequalities which were violated by x.
# """
function violated_facets(facets, x)
    violations = filter(
        facet -> facet_violation(facet, x) > 0 && !facet_contains(facet, x),
        facets
    )

    violations
end

# """
# facets_contain(facets, x)
#
#     Returns true if point x satisfies all facet inequalities.
#
# Inputs:
#     facet: Row Vector, contains hyperplane constant and normal vector of hyperplane.
#     x: Col Vector, a point the convex polytope space.
#
# Outputs:
#     x_is_contained: Boolean, true if x lies within or along the polytope boundary.
# """
function facets_contain(facets, x)
    x_is_contained = all(map(
        facet -> facet_violation(facet, x) <= 0,
        facets
    ))

    x_is_contained
end

# """
# facet_violation(facet, x)
#
#     Computes the violation quantity of point x for a given facet. A violation > 0
#     indicates a point which does not satisfy the facet inequality.
#
# Inputs:
#     facet: Row Vector, contains hyperplane constant and normal vector of hyperplane.
#     x: Col Vector, a point the convex polytope space.
#
# Outputs:
#     violation: Float, Inequality is satisfied for violation <= 0.
# """
function facet_violation(facet, x)
    violation = NaN
    if length(x) == length(facet)
        violation = (facet*x)[1]
    elseif (length(x) + 1) == length(facet)
        violation = facet[1] + (facet[2:end]'*x)[1]
    else
        throw(ArgumentError("dimension mismatch: facet not 'dual' to coordinate"))
    end

    violation
end

# """
# facet_distance(facet, x)
#
#     Computes the (shortest) distance between point x and a given facet. A
#     distance < 0 indicates that the point lies within the polytope.
#
# Inputs:
#     facet: Row Vector, contains hyperplane constant and normal vector of hyperplane.
#     x: Col Vector, a point the convex polytope space.
#
# Outputs:
#     dist: Float, distance a point is from a bounding hyperplane.
# """
function facet_distance(facet, x)
    dist = ConvexPolytope.facet_violation(facet, x)/norm(facet[2:end])
    dist
end

# """
# distance(x1, x2)
#     Computes the distance between two points in behavior space.
#
# Inputs:
#     x1/2: Col Vector, represents a point in behavior space.
#
# Output:
#     dist: Integer, the Euclidean distance between the two point.
# """
function vertex_distance(x1, x2)
    dist = norm(x1 - x2)
    dist
end

# """
# collect_vertices(facets, vertices):
#
#     For each facet, find all vertices that lie on that facet and construct a
#     dictionary linking the facet to the set of vertices.
#
# Input:
#     facets: Array, contains valid facet inequalities
#     vertices: Array, contains valid vertices
#
# Output:
#     Array, contains dictionaries in the following form
#         {
#             "facet": [facet inequality],
#             "vertices": [[v1],[v2],...,[vN]]
#         }
# """
function collect_vertices(facets, vertices)
    fv_groups = map(
        (f)->Dict(
            "facet" => f,
            "vertices" => filter(v->ConvexPolytope.facet_contains(f, v), vertices)
        ), facets)

    fv_groups
end

# """
# collect_facets(vertices, facets):
#
#     For each vertex, find all facets which that vertex lies on and construct
#     a dictionary linking each vertex to all of its facets.
#
# Input:
#     vertices: Array, contains valid vertices
#     facets: Array, contains valid facet inequalities
#
# Output:
#     Array, contains dictionaries of the following form
#         {
#             "vertex": [vertex]
#             "facets": [[f1],[f2],...,[fN]]
#         }
# """
function collect_facets(vertices, facets)
    vf_groups = map(
        (v)->Dict(
            "vertex" => v,
            "facets" => filter(f->ConvexPolytope.facet_contains(f, v), facets)
        ), vertices)

    vf_groups
end

end
