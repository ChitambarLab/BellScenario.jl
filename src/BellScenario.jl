module BellScenario

using QBase

# include external modules
include("./ConvexPolytope.jl")

# include internal modules
include("./LocalPolytope.jl")
include("./Behavior.jl")
include("./Degeneracy.jl")
include("./PrepareAndMeasure.jl")
include("./DichotomicLocalPolytope.jl")
include("./QuantumBehavior.jl")
include("./QuantumOpt.jl")

export QBase, ConvexPolytope, Degeneracy, LocalPolytope, Behavior, QuantumBehavior, QuantumOpt, PrepareAndMeasure

# function local_polytope()
#
# end
#
# function local_vertices()
# end
#
# function local_facets()
#
# end
#
# function optimize_quantum_nonlocality()
#
# end

end
