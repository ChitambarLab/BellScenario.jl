"""
*Compute the bounds of classical Bell scenarios.*

The `LocalPolytope.jl` module provides tools for characterizing the convex polyhedral
structure of the correlations attainable using classical resources in a Bell scenario.

Each distinct strategy for a given Bell scenario corresponds to a point in vector
space.
The complete set of attainable strategies for a particular Bell scenario is denoted
``\\mathbf{P}``.
Given the constraints of normalization and non-negativity, the extreme points of
``\\mathbf{P}`` correspond to deterministic strategies (matrices with 0,1 elements).
Furthermore, when shared randomness is used in a Bell scenario, any two strategies
can be mixed together in a convex combination.
Hence, ``\\mathbf{P}`` is a convex polyhedron referred to as the *local polytope*.

A convex polytope has two equivalent descriptions:

1. **V-Description**: The polytope is the convex hull of a set of extreme points
    known as vertices ``\\mathbf{V}``,
    ``\\quad\\mathbf{P} = \\text{conv}(\\mathbf{V})``.
2. **H-Description**: The polytope is the intersection of a set of linear half-spaces
    known as facets ``\\mathbf{F}``,
    ``\\quad\\mathbf{P} = \\cap\\mathbf{F}``.

These two descriptions are equivalent,

```math
\\mathbf{P} = \\text{conv}(\\mathbf{V}) = \\cap\\mathbf{F},
```

and there exists a transformation between the set of vertices and facets
``\\mathbf{V} \\leftrightarrow \\mathbf{F}``.

For a given Bell scenario, the V-Description is typically easier to compute,
however, the H-Description describes the bounds of the local polytope as
as linear inequalities known as Bell inequalities.
Bell inequalities are important because they provide a simple test to verify that
a strategy is not contained by the local polyotpe.
That is, if a Bell inequality is violated, the classical resources considered for
the local polytope are not sufficient to reproduce the strategy.
Hence, a Bell violation witnesses the use of resources with greater operational value.
Bell violations are often used to characterize the advantages of quantum resources,
however, it is important to note that a Bell violation can simply mean that more classical
resources were used than anticipated.

### Module Exports:
* [`vertices`](@ref) - Compute the set of extreme-points for the local polytope.
* [`num_vertices`](@ref) - The number of vertices for the local polytope.
* [`vrep`](@ref) - Construct a [`Polyhedron`](https://juliapolyhedra.github.io/Polyhedra.jl/stable/polyhedron/)
    in the vertex representation.
* [`facets`](@ref) - Compute the linear inequalities which bound the local polytope.
* [`generator_vertex`](@ref) - Provide a canonical form for a vertex.
* [`generator_facet`](@ref) - Provide a canonical form for a facet.
* [`adjacency_decomposition`](@ref) - Efficiently compute the generator facets for the local
    polytope using the adjacency decomposition technique.
"""
module LocalPolytope

using LinearAlgebra, Combinatorics
using XPORTA, Polyhedra

import Polyhedra: vrep

using ..BellScenario

include("./vertices.jl")

"""
    vrep(scenario::Scenario; vertices_kwargs...) :: XPORTA.Polyhedron

Constructs a `Polyhedron` using the vertex representation.
See [Polyhedra.jl](https://github.com/JuliaPolyhedra/Polyhedra.jl)
for more details.
The `vertices_kwargs` keyword arguments are passed through to the `vertices`
function for each `Scenario`.

!!! note "Return Type"
    This function differs from the Polyhedra.jl implementation in that it returns
    a `Polyhedron` type rather than a `V-Representation`.
    This is done to reduce the number of steps required to construct a polyhedron.
"""
function vrep(scenario::Scenario; vertices_kwargs...) :: XPORTA.Polyhedron
    polyhedron(vrep(vertices(scenario; vertices_kwargs...)), XPORTA.Library())
end

include("./facets.jl")
include("./generators.jl")
include("./adjacency_decomposition.jl")
include("./utils.jl")

# legacy code
include("../Legacy/LocalPolytope.jl")

end
