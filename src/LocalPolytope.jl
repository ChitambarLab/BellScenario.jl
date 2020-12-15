module LocalPolytope

# """
# Module Description:
#   This module houses functions used to compute local polytopes for bipartite
#   bell scenarios via the method of vertex (deterministic behavior) enumeration.
# """

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
