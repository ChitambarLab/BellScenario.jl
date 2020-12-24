"""
The `LocalPolytope` module characterizes the classical bounds of Bell scenarios.
The set of strategies attainable by classical devices forms a convex polyhedron
known as the *local polytope*.
Each distinct strategy represents a point in the local polytope and the extreme
points (vertices) of the local polytope are deterministic strategies.
When shared randomness is permitted between all black-boxes involved in the Bell
scenario, any mixture of deterministic strategies can be taken and the set forms
a convex polytope.
A convex polytope has two equivalent descriptions,

1. The V-Description: The polytope is the convex hull of a set of vertices, ``\\quad\\text{conv}(\\mathbf{V})``.
2. The H-Description: The polytope is the intersection of a set of linear half-spaces ``\\quad\\cap\\mathbf{F}``.

```math
\\text{conv}(\\mathbf{V}) = \\cap\\mathbf{F}
```

Typically, the V-Description is easier to compute, however, the H-Description
describes the bound as linear inequalities which provide a simple test for inclusion
or exclusion from the local polytope.

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
