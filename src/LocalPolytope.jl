"""
The `LocalPolytope` module is a submodule of [`BellScenario`](@ref). If shared
randomness is held between each of the black-boxes involved in the Bell scenario,
then the set of conditional probabilities attainable by the black-boxes form a
convex polytope regarded as the Local Polytope. A convex polytope has two equivalent
descriptions,

1. The V-Description: The polytope is the convex hull of a set of extreme-points.
2. The H-Description: The polytope is the intersection of a set of linear half-spaces.

The V-Description is typically easier to compute while the H-Description is useful
because it inherently provides a set of inequalities which confirm whether a point
is contained within the polytope. 

This module
exports the following memthods:
* [`vertices`](@ref) - Compute the set of extreme-points for the Local Polytope.
* [`facets`](@ref) - Compute the linear inequalities which bound the Local Polytope.
* [`generator_vertex`](@ref) - Provide a canonical form for a vertex.
* [`generator_facet`](@ref) - Provide a canonical form for a facet.
* [`adjacency_decomposition`](@ref) - Efficiently compute the generating facets for the Local
    Polytope using the adjacency decomposition technique.
"""
module LocalPolytope

using QBase: QMath
using LinearAlgebra
using XPORTA

using ..BellScenario

include("./LocalPolytope/vertices.jl")
include("./LocalPolytope/facets.jl")
include("./LocalPolytope/generators.jl")
include("./LocalPolytope/adjacency_decomposition.jl")
include("./LocalPolytope/utils.jl")

# legacy code
include("./Legacy/LocalPolytope.jl")

end
