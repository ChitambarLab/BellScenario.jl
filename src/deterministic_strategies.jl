# Types
export DeterministicStrategy

# Type validation methods
export is_deterministic

# constructor methods
export deterministic_strategies

"""
    DeterministicStrategy(conditionals :: Matrix{Int64}) <: AbstractStrategy{Int64}

A strategy matrix describing the behavior of a deterministic black-box.
A strategy *deterministic* if its elements satisfy ``P(y|x)\\in\\{0,1\\}`` in addition
to the non-negativity and normalization constraints.
By default, the constructor creates a strategy for a 'BlackBox' scenario, however,
a `Scenario` can be passed to the `DeterministicStrategy` constructor.

    DeterministicStrategy(conditionals :: Matrix{Int}, scenario :: Scenario)

The product of two deterministic strategies is a `DeterministicStrategy`.

    *(S1::DeterministicStrategy, S2::DeterministicStrategy) :: DeterministicStrategy

When multiplied a new `Scenario` can be specified.

    *(
        S1::DeterministicStrategy,
        S2::DeterministicStrategy,
        scenario::Scenario
    ) :: Deterministic Strategy

### Errors:

A `DomainError` is thrown if:
* The provided `Scenario` does not match the dimension of the `conditionals` matrix.
* The elements of `conditionals` are not `0` or `1`.
* The strategy elements are not non-negative and normalized.
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
    ) = is_deterministic(conditionals) ? new(conditionals, BlackBox(size(conditionals)...)) : throw(
        DomainError(conditionals, "conditionals are not a deterministic (`{0,1}`) distribution.")
    )
end

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
    is_deterministic( strategy :: AbstractMatrix  ) :: Bool

Returns `true` if all elements of `strategy` are either `0` or `1` and the matrix
is a valid conditional probability distribution.
"""
function is_deterministic(S :: AbstractMatrix) :: Bool
    is_conditionals = is_conditional_distribution(S)
    is_01 = all(i -> (S[i] == 0) || (S[i] == 1), 1:length(S))

    (is_conditionals && is_01)
end

"""
    deterministic_strategies(scenario :: BlackBox) :: Vector{Matrix{Int64}}

    deterministic_strategies(num_out :: Int64, num_in :: Int64) :: Vector{Matrix{Int64}}

Enumerates the set of deterministic strategies for the specified `BlackBox`. For
performance, enumerated deterministic strategies are left as matrices and not
constructed into [`DeterministicStrategy`](@ref) types.
"""
deterministic_strategies(
    num_out :: Int64,
    num_in :: Int64
) = deterministic_strategies(BlackBox(num_out,num_in))

function deterministic_strategies(scenario :: BlackBox) :: Vector{Matrix{Int64}}
    num_in = scenario.num_in
    num_out = scenario.num_out

    strategies = (num_out == 1) ? Matrix{Int64}[fill(1,(1,num_in))] : map( i -> begin
        m = zeros(Int64, num_out, num_in)
        base_n_array = digits(i, base = num_out, pad = num_in) .+ 1

        for j in 1:num_in
            m[base_n_array[j],j] = 1
        end

        m
    end, 0:(num_out^num_in - 1))

    strategies
end
