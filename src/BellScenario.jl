module BellScenario

using QBase, LinearAlgebra

import Base: *, convert

# exported modules
export LocalPolytope

# Legacy modules
export ConvexPolytope, Degeneracy, Behavior, QuantumBehavior, QuantumOpt

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

include("./LocalPolytope.jl")
using .LocalPolytope

# include internal modules
include("./Legacy/ConvexPolytope.jl")
using .ConvexPolytope

include("./Legacy/Behavior.jl")
using .Behavior
include("./Legacy/Degeneracy.jl")
using .Degeneracy

# include internal modules

include("./Legacy/DichotomicLocalPolytope.jl")
using .DichotomicLocalPolytope
#
include("./Legacy/QuantumBehavior.jl")
include("./Legacy/QuantumOpt.jl")
include("./Legacy/LocalSignaling.jl")

end
