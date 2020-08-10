export vertices, num_vertices

"""
    vertices(
        PM :: PrepareAndMeasure;
        rank_d_only = false :: Bool
    ) ::  Vector{Matrix{Int64}}

Generates the deterministic strategies for the local polytope of `PrepareAndMeasure`
scenarios.  The `rank_d_only` keyword arg specifies  whether to  exclude vertices
which use fewer dits  of communication and thus have a  matrix rank less than d.
"""
function vertices(PM :: PrepareAndMeasure; rank_d_only = false :: Bool) ::  Vector{Matrix{Int64}}
    vertices = Vector{Matrix{Int64}}(undef,0)
    lower_dits_bound = rank_d_only ? PM.d : 1
    for dits in lower_dits_bound:PM.d
        d_perms = QMath.permutations(1:dits)
        X_partitions = QMath.stirling2_partitions(PM.X, dits)
        B_combinations = QMath.combinations(1:PM.B, dits)

        for B in B_combinations
            for d in d_perms
                for X in X_partitions
                    # construct matrix m from permutation ids. This is faster
                    # than performing matrix multiplication
                    m = zeros(Int64, PM.B, PM.X)
                    for i in 1:dits
                        row_id = B[i]
                        partition_ids = X[d[i]]

                        m[row_id, partition_ids] .= 1
                    end
                    push!(vertices, m)
                end
            end
        end
    end

    vertices
end

"""
    num_vertices( PM :: PrepareAndMeasure ) :: Int64

Counts the numbr of polytope vertices for the specified `PrepareAndMeasure` scenario.
"""
function num_vertices(PM :: PrepareAndMeasure) :: Int64
    sum(map(i -> QMath.stirling2(PM.X, i)*binomial(PM.B, i)*factorial(i), 1:PM.d))
end
