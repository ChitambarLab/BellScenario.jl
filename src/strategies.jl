export AbstractStrategy, Strategy, DeterministicStrategy, strategy_dims, is_deterministic

"""
A stochastic matrix which represents a map from input to output for a given Bell scenario.
"""
abstract type AbstractStrategy{T} <: AbstractMatrix{T} end
Base.size(S::AbstractStrategy) = size(S.conditionals)
Base.getindex(S::AbstractStrategy, I::Vararg{Int,2}) = getindex(S.conditionals, I...)
Base.setindex!(S::AbstractStrategy, v, I::Vararg{Int,2}) = (S.conditionals[I...] = v)

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
        conditionals :: Matrix{<:Real},
        scenario :: Scenario
    ) =  size(conditionals) == strategy_dims(scenario) ? new(QMath.Conditionals(conditionals), scenario) : throw(
        DomainError(conditionals, "conditionals are not the correct dimension for the specified bell scenario.")
    )
    Strategy( conditionals :: Matrix{<:Real} ) = new(
        QMath.Conditionals(conditionals),
        BlackBox(reverse(size(conditionals))...)
    )
end

"""
A strategy matrix describing the deterministic behavior of a black-box.

    DeterministicStrategy(conditionals :: Matrix{Int64}) <: AbstractMatrix{Int64}

By default, the constructor creates a strategy for a 'BlackBox' scenario. However,
a `Scenario` can be passed to the `DeterministicStrategy` constructor.

    DeterministicStrategy(conditionals :: Matrix{Int}, scenario :: Scenario)

A `DomainError` is thrown if the provided `Scenario` does not match the dimension
of the `conditionals` matrix.

A `DomainError` is thrown if the elements of `conditionals` are not `0` or `1`.
"""
struct DeterministicStrategy <: AbstractStrategy{Int64}
    conditionals :: Matrix{Int64}
    scenario :: Scenario
    DeterministicStrategy(
        conditionals :: Matrix{<:Real},
        scenario :: Scenario
    ) = (is_deterministic(conditionals) && (size(conditionals) == strategy_dims(scenario))) ? new(conditionals, scenario) : throw(
        DomainError(conditionals, "conditionals are not a deterministic (`{0,1}`) distribution.")
    )
    DeterministicStrategy(
        conditionals :: Matrix{<:Real}
    ) = is_deterministic(conditionals) ? new(conditionals, BlackBox(reverse(size(conditionals))...)) : throw(
        DomainError(conditionals, "conditionals are not a deterministic (`{0,1}`) distribution.")
    )
end

"""
    is_deterministic( strategy :: AbstractMatrix  ) :: Bool

Returns `true` if all elements of `strategy` are either `0` or `1` and the matrix
is a valid conditional probability distribution.
"""
function is_deterministic(S :: AbstractMatrix) :: Bool
    is_conditionals = QMath.is_conditional_distribution(S)
    is_01 = all(i -> (S[i] == 0) || (S[i] == 1), 1:length(S))

    (is_conditionals && is_01)
end

"""
Two strategies can be multiplied together resulting in a new strategy, e.g. `S1*S2 = S3`.

    *(S1::AbstractStrategy, S2::AbstractStrategy) :: Strategy

The chained strategies can describe a new black-box scenario, therefore, `scenario`
parameter can also be passed to the multiplication operator.

    *(S1::AbstractStrategy, S2::AbstractStrategy, scenario::Scenario) :: Strategy

The product of two deterministic strategies is a `DeterministicStrategy`.

    *(S1::DeterministicStrategy, S2::DeterministicStrategy) :: DeterministicStrategy
"""
*(S1::AbstractStrategy, S2::AbstractStrategy) = Strategy(S1.conditionals*S2.conditionals)
*(
    S1::AbstractStrategy,
    S2::AbstractStrategy,
    scenario::Scenario
) = Strategy(S1.conditionals*S2.conditionals, scenario)

*(
    S1::DeterministicStrategy,
    S2::DeterministicStrategy
) = DeterministicStrategy(S1.conditionals*S2.conditionals)
*(
    S1::DeterministicStrategy,
    S2::DeterministicStrategy,
    scenario::Scenario
) = DeterministicStrategy(S1.conditionals*S2.conditionals, scenario)



"""
    strategy_dims( scenario :: Scenario ) :: Tuple{Int, Int}

Returns the dimensions of the `Conditionals` matrix describing a `Strategy` for
the `Scenario` at hand. Each `Scenario`, can place unique constraints on the matrix
dimensions, therfore, separate methods are called for each concrete `Scenario`.
"""
function strategy_dims(scenario::BlackBox) :: Tuple{Int, Int}
    (scenario.num_out, scenario.num_in)
end

function strategy_dims(scenario::PrepareAndMeasure) :: Tuple{Int, Int}
    (scenario.B, scenario.X)
end

function strategy_dims(scenario::Bipartite) :: Tuple{Int, Int}
    A = scenario.A
    B = scenario.B

    ( A.num_out * B.num_out, A.num_in * B.num_in )
end
