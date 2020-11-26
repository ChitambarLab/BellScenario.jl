export vertices, num_vertices

"""
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
    num_vertices( scenario :: LocalSignaling; rank_d_only = false :: Bool ) :: Int64

Counts the numbr of polytope vertices for the specified `LocalSignaling` scenario.
If `rank_d_only = true`, then only strategies using  `d`-dits are counted.
"""
function num_vertices(scenario :: LocalSignaling; rank_d_only = false :: Bool) :: Int64
    lower_dits_bound = rank_d_only ? scenario.d : 1
    sum(map(i -> QMath.stirling2(scenario.X, i)*binomial(scenario.B, i)*factorial(i), lower_dits_bound:scenario.d))
end
