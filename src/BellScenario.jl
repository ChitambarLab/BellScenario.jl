module BellScenario

using QBase

using Polyhedra: HalfSpace

import Base: *

export ConvexPolytope, Degeneracy, LocalPolytope, Behavior, Symmetry #, QuantumBehavior, QuantumOpt, PrepareAndMeasure

export rotate_facet, adjacent_facets, adjacency_decomposition

include("./scenarios.jl")
include("./strategies.jl")
include("./games.jl")

include("./file_io.jl")

"""
    facet_to_bell_game(num_in, num_out, gen_facet)

transforms the vector form of a facet to the minimal matrix with non-negative
elements.
"""
function facet_to_bell_game(gen_facet, scenario::Scenario) :: BellGame

    bound = -1*gen_facet[1]

    (num_out, num_in) = strategy_dims(scenario)

    matrix = behavior_to_strategy(num_in, num_out, gen_facet[2:end])

    new_matrix = zeros(num_out,num_in)
    for col_id in 1:num_in
        col = matrix[:,col_id]

        col_min = min(col...)

        new_col = col
        if col_min != 0
            new_col = -1*col_min .+ col

            bound = bound + -1*col_min
        end

        new_matrix[:,col_id] = new_col
    end

    (bound, new_matrix)
end

# """
# behavior_to_strategy(num_in, num_out, gen_behavior):
#
#     converts a behavior in the generalized representation into it's correpsonding
#     strategy matrix.
#
# Inputs:
#     num_in/out: Integer, bell scenario parameters (dimensions of strategies)
#     gen_behavior: Col Vector, a valid behavior in the generalized representation
#
# Output:
#     strategy: Matrix, the strategy isomorphic to the provided behavior.
# """
function behavior_to_strategy(gen_behavior, scenario::Scenario) :: Strategy
    behavior = []

    (num_out, num_in) = strategy_dims(scenario)

    if length(gen_behavior) == num_in * num_out
        behavior = gen_behavior
    elseif length(gen_behavior) - 1 == num_in * num_out
        behavior = gen_behavior[2:end]
    else
        throw(DomainError(gen_behavior, "./src/BellComm/LocalPolytope.jl: generalized behavior does not have proper dimensions"))
    end

    strategy = reshape(behavior, num_in, num_out)'

    strategy
end

# include internal modules
include("./ConvexPolytope.jl")
using .ConvexPolytope

include("./Behavior.jl")
using .Behavior
include("./Degeneracy.jl")
using .Degeneracy
include("./Symmetry.jl")
using .Symmetry
include("./PrepareAndMeasure.jl")
# include internal modules
include("./LocalPolytope.jl")
using .LocalPolytope
include("./DichotomicLocalPolytope.jl")
using .DichotomicLocalPolytope
#
include("./QuantumBehavior.jl")
include("./QuantumOpt.jl")

end
