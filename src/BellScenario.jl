"""
The `BellScenario` module provides:
1. A [`Scenario`](@ref) abstract type which serves as the parent of different
    multipartite black-box signaling scenarios.
2. An [`AbstractStrategy`](@ref) type to describe the contain the conditional probabilities
    of a given Bell Scenario.
3. An [`AbstractGame`](@ref) type to describe inequalities which bound the black-box statistics.
4. A quantum optimization library for tuning the parameters of quantum systems.
"""
module BellScenario

using QBase, LinearAlgebra

import Base: *, convert

# exported modules
export LocalPolytope, Nonlocality

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

# read/write and printing
include("./file_io.jl")

include("./LocalPolytope.jl")
using .LocalPolytope

include("./Nonlocality/Nonlocality.jl")
using .Nonlocality

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
