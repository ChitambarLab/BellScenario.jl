export vertices, num_vertices

export black_box_strategies

"""
    vertices(
        scenario :: BlackBox;
        rep = "normalized" :: String
    ) :: Vector{Vector{Int64}}

Generates the Local Polytope vertices for a BlackBox scenario. Valid represenations
are:

*`rep == "normalized"` or `rep == "generalized"`.

    vertices(
        scenario :: LocalSignaling;
        rep = "normalized" :: String
        rank_d_only = false :: Bool
    ) ::  Vector{Vector{Int64}}

Generates the deterministic strategies for the local polytope of `LocalSignaling`
scenarios.  The `rank_d_only` keyword arg specifies  whether to  exclude vertices
which use fewer dits  of communication and thus have a  matrix rank less than d.

!!! warning
    The vertices computed in this method are vectorized directly from a strategy
    matrix by column-majorization. These vertices are distinct from those produced
    by older `LocalPolytope.vertices()` methods which are row-majorized. This
    functionality is a trial run of performance improvements to polytope computation.
"""
function vertices(scenario :: LocalSignaling;
    rep = "normalized" :: String,
    rank_d_only = false :: Bool
) ::  Vector{Vector{Int64}}

    if !(rep in ("normalized", "generalized"))
        throw(DomainError(rep, "Argument `rep` must be either 'normalized' or 'generalized'."))
    end

    num_rows = (rep == "normalized") ? scenario.B - 1 : scenario.B
    lower_dits_bound = rank_d_only ? scenario.d : 1

    vertices = Vector{Vector{Int64}}(undef, num_vertices(scenario, rank_d_only = rank_d_only))
    vertex_id = 1
    for dits in lower_dits_bound:scenario.d
        d_perms = QMath.permutations(1:dits)
        X_partitions = QMath.stirling2_partitions(scenario.X, dits)
        B_combinations = QMath.combinations(1:scenario.B, dits)

        for B in B_combinations
            for d in d_perms
                for X in X_partitions
                    # construct matrix m from permutation ids. This is faster
                    # than performing matrix multiplication
                    m = zeros(Int64, num_rows, scenario.X)
                    for i in 1:dits
                        if B[i] > num_rows
                            continue
                        end

                        row_id = B[i]
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
see main docs block above for `vertices`
"""
function vertices(
    scenario :: BlackBox;
    rep = "normalized" :: String
) :: Vector{Vector{Int64}}
    if !(rep in ["normalized","generalized"])
        throw(DomainError(rep, "Argument `rep` must be either 'normalized' or 'generalized'."))
    end

    strategies = black_box_strategies(scenario)

    vertices = (rep == "normalized") ? map(
        s -> s[1:(end-1),:][:], strategies
    ) : map(
        s -> s[:], strategies
    )

    vertices
end

function vertices(scenario :: BipartiteNonSignaling)
end

"""
    num_vertices( scenario :: LocalSignaling; rank_d_only = false :: Bool ) :: Int64

Counts the numbr of polytope vertices for the specified `LocalSignaling` scenario.
If `rank_d_only = true`, then only strategies using  `d`-dits are counted.
"""
function num_vertices(scenario :: LocalSignaling; rank_d_only = false :: Bool) :: Int64
    lower_dits_bound = rank_d_only ? scenario.d : 1
    sum(map(i -> QMath.stirling2(scenario.X, i)*binomial(scenario.B, i)*factorial(i), lower_dits_bound:scenario.d))
end


"""
    black_box_strategies(scenario :: BlackBox) :: Vector{Matrix{Int64}}

    black_box_strategies(num_out :: Int64, num_in :: Int64) :: Vector{Matrix{Int64}}

Enumerates the set of deterministic strategies for the specified `BlackBox`.
"""
black_box_strategies(
    num_out :: Int64,
    num_in :: Int64
) = black_box_strategies(BlackBox(num_out,num_in))

function black_box_strategies(scenario :: BlackBox) :: Vector{Matrix{Int64}}
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
