module BellScenario

using QBase

import Base: *

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

export Scenario, BlackBox, Bipartite

export Strategy, Game

export strategy_dims

"""
An abstract type to represent general black-box scenarios.
"""
abstract type Scenario end

"""
    BlackBox(num_in :: Int64, num_out :: Int64) <: Scenario

A black-box scenario considering a single device. A black-box is an uncharacterized
device with a finite number of classical inputs and outputs.

A `DomainError` is throw if parameters `num_in` or `num_out` is less than 1.
"""
struct BlackBox <: Scenario
    num_in :: Int64
    num_out :: Int64
    BlackBox(
        num_in::Int64,
        num_out::Int64
    ) = ((num_in >= 1) && (num_out >= 1)) ? new(num_in, num_out) : throw(
        DomainError((num_in, num_out), "inputs must be ≥ 1")
    )
end

"""
    Bipartite(
        A :: Tuple{Int64, Int64},
        B :: Tuple{Int64, Int64};
        dits :: Int = 1,
        bidirectional :: Bool = false
    ) <: Scenario

A black-box scenario with two devices and possible communication between the
devices. The keyword parameter `dits` describes the number of dits used for
communication and `bidirectional` describes whether communication is simultaneous
between each party.
"""
struct Bipartite <: Scenario
    A :: BlackBox
    B :: BlackBox
    dits :: Int
    bidirectional :: Bool
    Bipartite(
        A::Tuple{Int,Int}, B::Tuple{Int,Int};
        dits::Int=1, bidirectional::Bool=false
    ) = (dits >= 1) ? new(BlackBox(A...), BlackBox(B...), dits, bidirectional) : throw(
            DomainError(dits, "communication `dits` must be ≥ 1.")
        )
end

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

"""
    Game(game::Matrix{T}, β::Real) <: AbstractMatrix{T}

A `Game` is played by a `Matrix` and `β is the bounding or winning score of the
game. A `Game` is matrix representation of a linear inequality. Each element is
a linear scale factor for an element of a strategy matrix.

Type parameter `T` is typically an `Int64` or `Float64`.
"""
struct Game{T} <: AbstractMatrix{T}
    game :: Matrix{T}
    Base.size(G::Game) = size(G.game)
    Base.getindex(G::Game, I::Vararg{Int64,2}) = getindex(G.game, I...)
    Base.setindex!(G::Game, v, I::Vararg{Int64,2}) = (G.game[I...] = v)

    β :: Real
    Game(game::Matrix{T}, β::Real) where T = new{T}(game, β)
end

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


end
