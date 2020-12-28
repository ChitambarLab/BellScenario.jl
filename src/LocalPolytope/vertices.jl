export vertices, num_vertices, vertex_dims

"""
    vertices( scenario :: BlackBox;
        rep = "normalized" :: String
    ) :: Vector{Vector{Int64}}

Generates the Local Polytope vertices for a BlackBox scenario. Valid represenations
are:

*`rep == "normalized"` or `rep == "generalized"`.
"""
function vertices(
    scenario :: BlackBox;
    rep = "normalized" :: String
) :: Vector{Vector{Int64}}
    if !(rep in ["normalized","generalized"])
        throw(DomainError(rep, "Argument `rep` must be either 'normalized' or 'generalized'."))
    end

    strategies = deterministic_strategies(scenario)

    vertices = (rep == "normalized") ? map(
        s -> s[1:(end-1),:][:], strategies
    ) : map(
        s -> s[:], strategies
    )

    vertices
end

"""
    vertices( scenario :: LocalSignaling;
        rep = "normalized" :: String
        rank_d_only = false :: Bool
    ) ::  Vector{Vector{Int64}}

Generates the deterministic strategies for the local polytope of `LocalSignaling`
scenarios.  The `rank_d_only` keyword arg specifies  whether to  exclude vertices
which use fewer dits  of communication and thus have a  matrix rank less than d.

!!! warning
    The vertices computed in this method are vectorized directly from a strategy
    matrix by column-majorization. These vertices are distinct from those produced
    by older `LocalPolytope.vertices()` methods which are row-majorized.
"""
function vertices(scenario :: LocalSignaling;
    rep = "normalized" :: String,
    rank_d_only = false :: Bool
) ::  Vector{Vector{Int64}}

    if !(rep in ("normalized", "generalized"))
        throw(DomainError(rep, "Argument `rep` must be either 'normalized' or 'generalized'."))
    end

    num_rows = (rep == "normalized") ? scenario.Y - 1 : scenario.Y
    lower_dits_bound = rank_d_only ? scenario.d : 1

    vertices = Vector{Vector{Int64}}(undef, num_vertices(scenario, rank_d_only = rank_d_only))
    vertex_id = 1
    for dits in lower_dits_bound:scenario.d
        d_perms = QMath.permutations(1:dits)
        X_partitions = QMath.stirling2_partitions(scenario.X, dits)
        Y_combinations = QMath.combinations(1:scenario.Y, dits)

        for Y in Y_combinations
            for d in d_perms
                for X in X_partitions
                    # construct matrix m from permutation ids. This is faster
                    # than performing matrix multiplication
                    m = zeros(Int64, num_rows, scenario.X)
                    for i in 1:dits
                        if Y[i] > num_rows
                            continue
                        end

                        row_id = Y[i]
                        partition_ids = X[d[i]]

                        m[row_id, partition_ids] .= 1
                    end
                    vertices[vertex_id] = m[:]

                    vertex_id +=  1
                end
            end
        end
    end

    vertices
end

"""
    vertices( scenario :: BipartiteNonSignaling,
        rep="non-signaling" :: String
    ) :: Vector{Vector{Int64}}

Enumerates the LocalPolytope vertices for the [`BipartiteNonSignaling`](@ref) scenario.
Valid representations for the vertices include:
* `"non-signaling"`, `"normalized"`, `"generalized"`

A `DomainError` is thrown if a valid representation is not specified.
"""
function vertices(scenario :: BipartiteNonSignaling, rep="non-signaling" :: String) :: Vector{Vector{Int64}}
    alice_strategies = deterministic_strategies(scenario.A, scenario.X)
    bob_strategies = deterministic_strategies(scenario.B, scenario.Y)

    if rep == "non-signaling"
        α_strategies = map(s -> s[1:end-1,:], alice_strategies)
        β_strategies = map(s -> s[1:end-1,:], bob_strategies)

        dim_α = scenario.X*(scenario.A - 1)
        dim_β = scenario.Y*(scenario.B - 1)
        dim_v = dim_α + dim_β + dim_α*dim_β

        strategies = map( (α, β) for α in α_strategies for β in β_strategies) do (α,β)
            v = zeros(Int64, dim_v)

            v[1:dim_α+dim_β] = vcat(α[:],β[:])
            v[dim_α+dim_β+1:end] = kron(α,β)[:]

            v
        end

        return strategies

    elseif rep in ("normalized", "generalized")

        strategies = (rep == "normalized") ? map(
            (α, β) for α in alice_strategies for β in bob_strategies) do (α,β)
                kron(α,β)[1:end-1,:][:]
            end : map(
            (α, β) for α in alice_strategies for β in bob_strategies) do (α,β)
                kron(α,β)[:]
            end

        return strategies
    else
        throw(DomainError(rep, "`rep in (\"non-signaling\",\"normalized\",\"generalized\")` must be satisfied" ))
    end
end

"""
    num_vertices( scenario :: BlackBox ) :: Int64

For ``n`` outputs and ``m`` inputs the number of vertices ``|\\mathcal{V}|`` are counted:

```math
|\\mathcal{V}| = n^m
```
"""
function num_vertices(scenario :: BlackBox) :: Int64
    scenario.num_out^scenario.num_in
end

"""
    num_vertices( scenario :: LocalSignaling;
        rank_d_only = false :: Bool
    ) :: Int64

If `rank_d_only = true`, then only strategies using  `d`-dits are counted. For
``X`` inputs and ``Y`` outputs the number of vertices ``|\\mathcal{V}|`` are counted:

```math
|\\mathcal{V}| = \\sum_{c=1}^d \\left\\{X \\atop c \\right\\}\\binom{Y}{c}c!
```
"""
function num_vertices(scenario :: LocalSignaling; rank_d_only = false :: Bool) :: Int64
    lower_dits_bound = rank_d_only ? scenario.d : 1
    sum(map(i -> QMath.stirling2(scenario.X, i)*binomial(scenario.Y, i)*factorial(i), lower_dits_bound:scenario.d))
end

"""
    num_vertices( scenario :: BipartiteNonSignaling ) :: Int64

For two non-signaling black-boxes with ``X`` and ``Y`` inputs and ``A`` and ``B``
outputs respectively, the number of vertices ``|\\mathcal{V}|`` are counted:

```math
|\\mathcal{V}| = A^X B^Y
```
"""
function num_vertices(scenario :: BipartiteNonSignaling) :: Int64
    black_box_A = BlackBox(scenario.A, scenario.X)
    black_box_B = BlackBox(scenario.B, scenario.Y)

    num_vertices(black_box_A)*num_vertices(black_box_B)
end

"""
For the given [`Scenario`](@ref), returns the length of the vertex in the representation
specified by `rep`. A `DomainError` is thrown if the `rep` is invalid.

    vertex_dims(scenario :: Union{BlackBox,LocalSignaling}, rep :: String) :: Int64

Valid values of `rep` are `"normalized"` and `"generalized"`.
"""
function vertex_dims(scenario :: Union{BlackBox,LocalSignaling}, rep :: String) :: Int64
    if !(rep in ["normalized", "generalized"])
        throw(DomainError(rep, "Argument `rep` must be either \"normalized\" or \"generalized\"."))
    end

    strat_dims = strategy_dims(scenario)

    (rep == "normalized") ? strat_dims[2] * (strat_dims[1] - 1) : strat_dims[2] * strat_dims[1]
end

"""
    vertex_dims( scenario:: BipartiteNonSignaling, rep :: String ) :: Int64

Valid values for `rep` include:
* "non-signaling"
* "normalized"
* "generalized"
"""
function vertex_dims(scenario :: BipartiteNonSignaling, rep :: String) :: Int64
    if !(rep in ["non-signaling", "normalized", "generalized"])
        throw(DomainError(rep, "Argument `rep` must be either \"non-signaling\", \"normalized\", or \"generalized\"."))
    end

    dim_α = scenario.X*(scenario.A - 1)
    dim_β = scenario.Y*(scenario.B - 1)

    (rep == "non-signaling") ? dim_α + dim_β + dim_α*dim_β : vertex_dims(
        BlackBox(strategy_dims(scenario)...), rep)
end
