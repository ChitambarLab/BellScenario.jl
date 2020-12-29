"""
*Types and constructors that represent Bell scenarios, their statistics, and their bounds.*

The BellScenario.jl module is the base library for Bell scenario analysis.

### Features:
* [`Scenario`](@ref): abstract type describing Bell scenarios.
* [`AbstractStrategy`](@ref): abstract type for Bell scenario statistics.
* [`AbstractGame`](@ref): abstract type for Bell scenario bounds.

A game theoretic framework is used to evaluate the performance of black-box systems.
A cost function testing the black-box system is regarded as a *game* while the statistics
generated by the black-box system are regarded as a *strategy*.
A strategy is played against a game to achieve a score, hence, this framework
provides a quantitative metric of the performance of a black-box system
in regards to a particular task.
"""
module BellScenario

using QBase, LinearAlgebra

import Base: *, convert

# exported modules
export LocalPolytope, Nonlocality

# Legacy modules
export ConvexPolytope, Degeneracy, Behavior, QuantumBehavior, QuantumOpt

# Scenarios
include("./scenarios.jl")

# Strategies
include("./strategies.jl")
include("./deterministic_strategies.jl")
include("./strategy_conversions.jl")
include("./quantum_strategies.jl")

# Games
include("./games.jl")
include("./game_conversions.jl")

# read/write and printing
include("./file_io.jl")

include("./LocalPolytope/LocalPolytope.jl")
using .LocalPolytope

include("./Nonlocality/Nonlocality.jl")
using .Nonlocality

# include legacy modules
include("./Legacy/ConvexPolytope.jl")
using .ConvexPolytope

include("./Legacy/Behavior.jl")
using .Behavior

include("./Legacy/Degeneracy.jl")
using .Degeneracy

include("./Legacy/DichotomicLocalPolytope.jl")
using .DichotomicLocalPolytope

include("./Legacy/QuantumBehavior.jl")
include("./Legacy/QuantumOpt.jl")
include("./Legacy/LocalSignaling.jl")

end
