module BellScenario

using QBase, LinearAlgebra

import Base: *, convert

export ConvexPolytope, Degeneracy, LocalPolytope, Behavior, Symmetry #, QuantumBehavior, QuantumOpt, LocalSignaling

# black-box scenarios
include("./scenarios.jl")

# fundamental data structures
include("./strategies.jl")
include("./strategy_utilities.jl")
include("./deterministic_strategies.jl")

include("./games.jl")
include("./game_conversions.jl")

# quantum scenarios
include("./quantum_strategies.jl")
include("./quantum_opt.jl")

# read/write and printing
include("./file_io.jl")

# include internal modules
include("./ConvexPolytope.jl")
using .ConvexPolytope

include("./Behavior.jl")
using .Behavior
include("./Degeneracy.jl")
using .Degeneracy
include("./Symmetry.jl")
using .Symmetry
include("./LocalSignaling.jl")
# include internal modules
include("./LocalPolytope.jl")
using .LocalPolytope
include("./DichotomicLocalPolytope.jl")
using .DichotomicLocalPolytope
#
include("./QuantumBehavior.jl")
include("./QuantumOpt.jl")

end
