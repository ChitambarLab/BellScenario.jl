module BellScenario

using QBase

using Polyhedra: HalfSpace

import Base: *

export ConvexPolytope, Degeneracy, LocalPolytope, Behavior, Symmetry #, QuantumBehavior, QuantumOpt, PrepareAndMeasure

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


export Strategy, AbstractGame, Game, BellGame

export strategy_dims

export rotate_facet, adjacent_facets, adjacency_decomposition
include("./scenarios.jl")

"""
A strategy matrix describing the statistical behavior of a black-box.

    Strategy(conditionals :: Matrix{<:Real}) <: AbstractMatrix{Float64}

By default, the constructor creates a strategy for a 'BlackBox' scenario. However,
a `Scenario` can also be passed to the `Strategy` constructor.

    Strategy(conditionals :: Matrix{<:Real}, scenario :: Scenario)

A `DomainError` is thrown if the provided `Scenario` does not match the dimension
of the `conditionals` matrix.
"""
struct Strategy <: AbstractMatrix{Float64}
    conditionals ::  QMath.Conditionals
    Base.size(S::Strategy) = size(S.conditionals)
    Base.getindex(S::Strategy, I::Vararg{Int,2}) = getindex(S.conditionals, I...)
    Base.setindex!(S::Strategy, v, I::Vararg{Int,2}) = (S.conditionals[I...] = v)

    scenario :: Scenario
    Strategy(
        conditionals :: Matrix{<:Real},
        scenario::Scenario
    ) =  size(conditionals) == strategy_dims(scenario) ? new(QMath.Conditionals(conditionals), scenario) : throw(
        DomainError(conditionals, "conditionals are not the correct dimension for the specified bell scenario.")
    )
    Strategy( conditionals :: Matrix{<:Real} ) = new(
        QMath.Conditionals(conditionals),
        BlackBox(reverse(size(conditionals))...)
    )
end

"""
Two strategies can be multiplied together resulting in a new strategy, e.g. `S1*S2 = S3`.

    *(S1::Strategy, S2::Strategy) :: Strategy

The chained strategies can describe a new black-box scenario, therefore, `scenario`
parameter can also be passed to the multiplication operator.

    *(S1::Strategy, S2::Strategy, scenario::Scenario) :: Strategy
"""
*(S1::Strategy, S2::Strategy) = Strategy(S1.conditionals*S2.conditionals)
*(S1::Strategy, S2::Strategy, scenario::Scenario) = Strategy(S1.conditionals*S2.conditionals, scenario)

"""
    strategy_dims( scenario :: Scenario ) :: Tuple{Int, Int}

Returns the dimensions of the `Conditionals` matrix describing a `Strategy` for
the `Scenario` at hand. Each `Scenario`, can place unique constraints on the matrix
dimensions, therfore, separate methods are called for each concrete `Scenario`.
"""
function strategy_dims(scenario::BlackBox) :: Tuple{Int, Int}
    (scenario.num_out, scenario.num_in)
end

function strategy_dims(scenario::Bipartite) :: Tuple{Int, Int}
    A = scenario.A
    B = scenario.B

    ( A.num_out * B.num_out, A.num_in * B.num_in )
end

abstract type AbstractGame{T} <: AbstractMatrix{T} end
Base.size(G::AbstractGame) = size(G.game)
Base.getindex(G::AbstractGame, I::Vararg{Int64,2}) = getindex(G.game, I...)
Base.setindex!(G::AbstractGame, v, I::Vararg{Int64,2}) = (G.game[I...] = v)

"""
    Game(game::Matrix{T}, β::Real) <: AbstractMatrix{T}

A `Game` is played by a `Matrix` and `β is the bounding or winning score of the
game. A `Game` is matrix representation of a linear inequality. Each element is
a linear scale factor for an element of a strategy matrix.

Type parameter `T` is typically an `Int64` or `Float64`.
"""
struct Game{T} <: AbstractGame{T}
    game :: Matrix{T}
    β :: Real
    scenario :: Scenario
    Game(game::Matrix{T}, β::Real) where T = new{T}(game, β)
end

struct BellGame <: AbstractGame{Int64}
    game :: Matrix{Int64}
    β :: Int64
    scenario :: Scenario
    BellGame(game::Matrix{Int64}, β::Int64) = new(game, β)
end

"""
    facet_to_matrix(num_in, num_out, gen_facet)

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
#include("./PrepareAndMeasure.jl")
# include internal modules
include("./LocalPolytope.jl")
using .LocalPolytope
include("./DichotomicLocalPolytope.jl")
using .DichotomicLocalPolytope
#
# include("./QuantumBehavior.jl")
# include("./QuantumOpt.jl")

end
