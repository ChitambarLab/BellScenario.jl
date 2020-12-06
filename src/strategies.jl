# Types
export AbstractStrategy, Strategy

# Validation methods
export strategy_dims

# Constructors
export random_strategy

using Random: shuffle!

"""
A stochastic matrix which represents a map from input to output for a given Bell scenario.

# Base library extensions for  `AbstractStrategy`:

Two strategies can be multiplied together resulting in a new strategy, e.g. `S1*S2 = S3`.

    *(S1::AbstractStrategy, S2::AbstractStrategy) :: Strategy

The chained strategies can describe a new black-box scenario, therefore, `scenario`
parameter can also be passed to the multiplication operator.

    *(S1::AbstractStrategy, S2::AbstractStrategy, scenario::Scenario) :: Strategy
"""
abstract type AbstractStrategy{T} <: AbstractMatrix{T} end
Base.size(S::AbstractStrategy) = size(S.conditionals)
Base.getindex(S::AbstractStrategy, I::Vararg{Int,2}) = getindex(S.conditionals, I...)
Base.setindex!(S::AbstractStrategy, v, I::Vararg{Int,2}) = (S.conditionals[I...] = v)

*(S1::AbstractStrategy, S2::AbstractStrategy) = Strategy(S1.conditionals*S2.conditionals)
*(
    S1::AbstractStrategy,
    S2::AbstractStrategy,
    scenario::Scenario
) = Strategy(S1.conditionals*S2.conditionals, scenario)

"""
A strategy matrix describing the statistical behavior of a black-box.

    Strategy(conditionals :: Matrix{<:Real}) <: AbstractMatrix{Float64}

By default, the constructor creates a strategy for a 'BlackBox' scenario. However,
a `Scenario` can be passed to the `Strategy` constructor.

    Strategy(conditionals :: Matrix{<:Real}, scenario :: Scenario)

A `DomainError` is thrown if the provided `Scenario` does not match the dimension
of the `conditionals` matrix.
"""
struct Strategy <: AbstractStrategy{Float64}
    conditionals ::  QMath.Conditionals
    scenario :: Scenario
    Strategy(
        conditionals :: QMath.Conditionals,
        scenario :: Scenario
    ) =  size(conditionals) == strategy_dims(scenario) ? new(conditionals, scenario) : throw(
        DomainError(conditionals, "conditionals are not the correct dimension for the specified bell scenario.")
    )
    Strategy(
        conditionals :: Matrix{<:Real},
        scenario :: Scenario
    ) = Strategy(QMath.Conditionals(conditionals), scenario)
    Strategy( conditionals :: Matrix{<:Real} ) = new(
        QMath.Conditionals(conditionals),
        BlackBox(size(conditionals)...)
    )
end

"""
    random_strategy(
        num_inputs :: Int64,
        num_outputs :: Int64;
    ) :: Strategy

Constructs a randomized [`Strategy`](@ref) matrix. Zeros are inserted into the
strategy to ensure that it does not closely resemble a uniform distribution.
"""
function random_strategy(num_inputs :: Int64, num_outputs :: Int64) :: Strategy
    prob_vecs = map(i -> begin
        num_nonzero_elements = rand(1:num_outputs)
        v = zeros(Float64, num_outputs)
        v[1:num_nonzero_elements] = rand(Float64, num_nonzero_elements)
        shuffle!(v)/norm(v,1)
    end, 1:num_inputs)

    Strategy(hcat(prob_vecs...))
end

"""
    strategy_dims( scenario :: Scenario ) :: Tuple{Int, Int}

Returns the dimensions of the `Conditionals` matrix describing a `Strategy` for
the `Scenario` at hand. Each `Scenario`, can place unique constraints on the matrix
dimensions, therfore, separate methods are called for each concrete `Scenario`.
"""
function strategy_dims(scenario::BlackBox) :: Tuple{Int64, Int64}
    (scenario.num_out, scenario.num_in)
end

function strategy_dims(scenario::LocalSignaling) :: Tuple{Int64, Int64}
    (scenario.B, scenario.X)
end

function strategy_dims(scenario::BipartiteNoSignaling) :: Tuple{Int64,Int64}
    (scenario.A * scenario.B, scenario.X * scenario.Y)
end

function strategy_dims(scenario::Bipartite) :: Tuple{Int, Int}
    A = scenario.A
    B = scenario.B

    ( A.num_out * B.num_out, A.num_in * B.num_in )
end
