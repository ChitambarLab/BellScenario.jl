# Types
export AbstractStrategy, Strategy

# Validation methods
export strategy_dims

# Constructors
export random_strategy

using Random: shuffle!

"""
An `AbstractStrategy` is an abstract type parent to all strategy matrices.
A black-box device is characterized by its conditional probabilities which can be
organized into a strategy matrix or column-stochastic map ``S : \\mathcal{X}
\\to \\mathcal{Y}``.
The strategy matrix ``S`` is constructed as follows,

```math
S = \\sum_{x,y} P(y|x) |y\\rangle\\langle x|
```

where the elements of a strategy matrix must be non-negative and normalized:
* ``P(y|x) \\geq 0``
* ``\\sum_y P(y|x) = 1``

Since strategies are just stochastic matrices, the product of two strategies yields
a new strategy, e.g. ``S_A*S_B = S_C``.

    *(S1::AbstractStrategy, S2::AbstractStrategy) :: Strategy

When multiplied together, two strategies may represent a new `Scenario`, hence,
a the multiplication operator can also be passed a `Scenario`.

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
    Strategy(conditionals :: Matrix{<:Real}) <: AbstractMatrix{Float64}

    Strategy(conditionals :: QMath.Conditionals) <: AbstractMatrix{Float64}

The `conditionals` parameter is a column stochastic matrix which can be provided
in a raw `Matrix{<:Real}` format or as a `QMath.Conditionals` type from the QBase.jl
package.
By default, the constructor creates a strategy for a 'BlackBox' scenario. However,
a `Scenario` can be passed to the `Strategy` constructor.

    Strategy(conditionals :: Matrix{<:Real}, scenario :: Scenario)

### Errors:

A `DomainError` is thrown if:
* The provided `Scenario` does not match the expected dimension of the `conditionals` matrix.
* The `conditionals` matrix is not a valid stochastic matrix, e.g. non-negative and normalized.
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

Constructs a randomized [`Strategy`](@ref) matrix.
Zeros are inserted into the strategy to ensure that it does not closely resemble
a uniform distribution.
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
    strategy_dims( scenario :: Scenario ) :: Tuple{Int64, Int64}

Returns the dimensions of the `QMath.Conditionals` matrix describing a `Strategy` for
the `Scenario` at hand. Each `Scenario`, can place unique constraints on the matrix
dimensions, therfore, separate methods are called for each concrete `Scenario`.
"""
function strategy_dims(scenario::BlackBox) :: Tuple{Int64, Int64}
    (scenario.num_out, scenario.num_in)
end

function strategy_dims(scenario::LocalSignaling) :: Tuple{Int64, Int64}
    (scenario.Y, scenario.X)
end

function strategy_dims(scenario::BipartiteNonSignaling) :: Tuple{Int64, Int64}
    (scenario.A * scenario.B, scenario.X * scenario.Y)
end

function strategy_dims(scenario::BipartiteSignaling) :: Tuple{Int64, Int64}
    A = scenario.A
    B = scenario.B

    ( A.num_out * B.num_out, A.num_in * B.num_in )
end
