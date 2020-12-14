# Types
export DeterministicStrategy

# Type validation methods
export is_deterministic

# constructor methods
export deterministic_strategies

"""
A strategy matrix describing the deterministic behavior of a black-box.

    DeterministicStrategy(conditionals :: Matrix{Int64}) <: AbstractMatrix{Int64}

By default, the constructor creates a strategy for a 'BlackBox' scenario. However,
a `Scenario` can be passed to the `DeterministicStrategy` constructor.

    DeterministicStrategy(conditionals :: Matrix{Int}, scenario :: Scenario)

A `DomainError` is thrown if the provided `Scenario` does not match the dimension
of the `conditionals` matrix.

A `DomainError` is thrown if the elements of `conditionals` are not `0` or `1`.

# Base library Extensions for `DeterministicStrategy`:

The product of two deterministic strategies is a `DeterministicStrategy`.

    *(S1::DeterministicStrategy, S2::DeterministicStrategy) :: DeterministicStrategy

Deterministic strategies correspond to vertices of the local polytope. A vertex is
simply represented by a `Vector{Int64}` and the `rep` argument specifies whether
the vertex is in the `"normalized"` or `"generalized"` representation.

Vertex (`Vector{Int64}`) -> `DeterministicStrategy`

    convert(
        ::Type{DeterministicStrategy},
        vertex  :: Vector{Int64},
        scenario :: Scenario;
        rep = "normalized" :: String
    )

`DeterministicStrategy` -> Vertex (`Vector{Int64}`)

    convert(
        ::Type{Vector{Int64}},
        strategy :: DeterministicStrategy;
        rep = "normalized" :: String
    )
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

function convert(::Type{DeterministicStrategy}, vertex::Vector{Int64}, scenario::Scenario; rep="normalized" :: String)
    if !(rep in ("normalized", "generalized"))
        throw(DomainError(rep, "Argument `rep` must be either 'normalized' or 'generalized'"))
    end

    s_dims = strategy_dims(scenario)

    strategy_matrix = (rep == "normalized") ? cat(
        reshape(vertex, (s_dims[1]-1, s_dims[2])),
        zeros(Int64, (1,s_dims[2])),
        dims = 1
    ) : reshape(vertex, s_dims)

    # for vertices in the normalized representation, if a column does not contain
    # a one, add it to the last row.
    if rep == "normalized"
        for col_id in 1:s_dims[2]
            if sum(strategy_matrix[:,col_id]) == 0
                strategy_matrix[s_dims[1], col_id] = 1
            end
        end
    end

    DeterministicStrategy(strategy_matrix, scenario)
end

function convert(::Type{Vector{Int64}}, strategy::DeterministicStrategy; rep = "normalized"::String)
    if !(rep in ("normalized", "generalized"))
        throw(DomainError(rep, "Argument `rep` must be either 'normalized' or 'generalized'"))
    end

    max_row = (rep == "normalized") ? size(strategy,1) - 1 : size(strategy,1)

    strategy[1:max_row,:][:]
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