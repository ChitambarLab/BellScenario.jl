export stirling2, stirling2_partitions, stirling2_matrices
export permutation_matrices, n_choose_k_matrices
export base_n_val

"""
    stirling2( n :: Int64, k :: Int64  ) :: Int64

Counts the number of ways to partition `n` items into `k` unlabelled groups.
This quantity is known as Stirling's number of the 2nd kind:

```math
\\left\\{n \\atop k \\right\\} = \\frac{1}{k!}\\sum_{i=0}^k (-1)^i\\binom{k}{i}(k-i)^n
```

Throws a `DomainError` if inputs do not satisfy `n ≥ k ≥ 1`.
"""
function stirling2(n :: Int64, k :: Int64) :: Int64
    if !(n ≥ k ≥ 1)
        throw(DomainError((n,k), "Inputs (n,k) do not satisfy `n ≥ k ≥ 1`."))
    end

    sum(map( i -> (-1)^i*binomial(k, i)*(k-i)^n, 0:k) )/factorial(k)
end

"""
    stirling2_partitions( n :: Int64, k :: Int64 ) :: Vector{Vector{Vector{Int64}}}

Enumerates the unique partitions of `n` items into `k` unlabelled sets.
Each partition is a vector containing a  set of `k` vectors designating each group.

E.g.

```jldoctest
julia> stirling2_partitions(4, 2)
7-element Vector{Vector{Vector{Int64}}}:
 [[1, 2, 3], [4]]
 [[3], [1, 2, 4]]
 [[1, 2], [3, 4]]
 [[1, 3], [2, 4]]
 [[2], [1, 3, 4]]
 [[2, 3], [1, 4]]
 [[1], [2, 3, 4]]
```

This recursive algorithm was inspired by [this blog](https://devblogs.microsoft.com/oldnewthing/20140324-00/?p=1413).
"""
function stirling2_partitions(n :: Int64, k :: Int64) :: Vector{Vector{Vector{Int64}}}
    if (n ==  0) || (k == 0) || (n < k)
        return Vector{Vector{Vector{Int64}}}(undef,0)
    elseif (k == 1)
        return [[[1:n...]]]
    else
        # partitions include all (n-1, k-1) paritions with [n] as k-th group
        partitions1 = _stirling2_add_partition(n, stirling2_partitions(n-1, k-1))

        # partitions include all (n-1,  k) partitions with [n] appended to each existing group
        partitions2 = _stirling2_extend_partitions(n, k, stirling2_partitions(n-1, k))

        return cat(partitions1, partitions2, dims=1)
    end
end

# helper function for adding a new group to a partition
function _stirling2_add_partition(n :: Int64, partitions :: Vector{Vector{Vector{Int64}}}) :: Vector{Vector{Vector{Int64}}}
    if length(partitions) == 0
        return [[[n]]] :: Vector{Vector{Vector{Int64}}}
    else
        return map(s -> cat(s, [[n]], dims=1), partitions) :: Vector{Vector{Vector{Int64}}}
    end
end

# helper function for extending existing groups with a new element
function _stirling2_extend_partitions(n :: Int64, k :: Int64, partitions :: Vector{Vector{Vector{Int64}}}) :: Vector{Vector{Vector{Int64}}}
    new_partitions = Vector{Vector{Vector{Int64}}}(undef, length(partitions)*k)
    id = 1

    _extend_partition_set(partition :: Vector{Vector{Int64}}) = begin
        for i in 1:k
            new_partitions[id] = cat(partition[1:k .!= i], [cat(partition[i], [n], dims=1)], dims=1)
            id += 1
        end
    end

    map(_extend_partition_set, partitions)

    new_partitions
end

"""
    stirling2_matrices( n :: Int64, k :: Int64 ) :: Vector{Matrix{Bool}}

Generates the set of matrices with `k` rows and `n` columns where rows correspond
to the groups and columns are the grouped elements. A non-zero element designates
that the column id is grouped into the corresponding row.

E.g.

```jldoctest
julia> stirling2_matrices(4, 2)
7-element Vector{Matrix{Bool}}:
 [1 1 1 0; 0 0 0 1]
 [0 0 1 0; 1 1 0 1]
 [1 1 0 0; 0 0 1 1]
 [1 0 1 0; 0 1 0 1]
 [0 1 0 0; 1 0 1 1]
 [0 1 1 0; 1 0 0 1]
 [1 0 0 0; 0 1 1 1]
```

A `DomainError` is thrown if `n ≥ k ≥ 1` is not satisfied.
"""
function stirling2_matrices(n :: Int64, k :: Int64) :: Vector{Matrix{Bool}}
    if !(n ≥ k ≥ 1)
        throw(DomainError((n,k), "Inputs (n,k) do not satisfy `n ≥ k ≥ 1`."))
    end

    _partition_to_matrix(partition :: Vector{Vector{Int64}}) :: Matrix{Bool} = begin
        m = fill(false, k, n)
        for row_id in 1:k
            col_ids = partition[row_id]
            m[row_id,col_ids] .= true
        end

        return m
    end

    map(_partition_to_matrix, stirling2_partitions(n, k))
end

"""
    permutation_matrices( dim :: Int64 ) :: Vector{Matrix{Bool}}

Generates the set of square permutation matrices of dimension `dim`.

E.g.

```jldoctest
julia> permutation_matrices(3)
6-element Vector{Matrix{Bool}}:
 [1 0 0; 0 1 0; 0 0 1]
 [1 0 0; 0 0 1; 0 1 0]
 [0 1 0; 1 0 0; 0 0 1]
 [0 0 1; 1 0 0; 0 1 0]
 [0 1 0; 0 0 1; 1 0 0]
 [0 0 1; 0 1 0; 1 0 0]
```
"""
function permutation_matrices(dim :: Int64) :: Vector{Matrix{Bool}}
    map(perm_ids -> Matrix{Bool}(I,dim,dim)[:,perm_ids], permutations(1:dim))
end

"""
    n_choose_k_matrices( n :: Int64, k :: Int64 ) :: Vector{Matrix{Bool}}

Generates a set of `n` by `k` matrices which represent all combinations of selecting
`k` columns from `n` rows.  Each column, contains a single non-zero element and
`k` rows contain a non-zero element.

E.g.

```jldoctest
julia> n_choose_k_matrices( 4, 2 )
6-element Vector{Matrix{Bool}}:
 [1 0; 0 1; 0 0; 0 0]
 [1 0; 0 0; 0 1; 0 0]
 [1 0; 0 0; 0 0; 0 1]
 [0 0; 1 0; 0 1; 0 0]
 [0 0; 1 0; 0 0; 0 1]
 [0 0; 0 0; 1 0; 0 1]
```

A `DomainError` is thrown if `n ≥ k ≥ 1` is not satisfied.
"""
function n_choose_k_matrices(n :: Int64, k :: Int64) :: Vector{Matrix{Bool}}
    if !(n ≥ k ≥ 1)
        throw(DomainError((n,k), "Inputs (n,k) must satisfy `n ≥ k ≥ 1`."))
    end
    _combo_to_matrix(combo :: Vector{Int64}) :: Matrix{Bool} = begin
        m = fill(false,n,k)
        for col_id in 1:k
            m[combo[col_id],col_id] = true
        end

        return m
    end

    map(_combo_to_matrix, combinations(1:n, k))
end

"""
    base_n_val(
        num_array :: Vector{Int64},
        base :: Int64;
        big_endian=true :: Bool
    ) :: Int64

Given an array representing a number in base-n returns the value of that
number in base-10.

Inputs:
* `num_array` - Vector containing semi-positive integers less than base.
* `base` - The base-n number represented by num_array.
* `big_endian` - `true` if most significant place is at index 1, else `false`.
"""
function base_n_val(num_array :: Vector{Int64}, base :: Int64; big_endian = true :: Bool) :: Int64
    if maximum(num_array) >= base
        throw(ArgumentError("max value of $num_array is not valid in $base"))
    elseif minimum(num_array) < 0
        throw(ArgumentError("min value of $num_array must be zero or greater"))
    end

    pad = length(num_array)

    place_vals = big_endian ? map(
        x -> base^(pad - x[1]), enumerate(num_array)
    ) : map(
        x -> base^(x[1] - 1), enumerate(num_array)
    )

    val = place_vals'*num_array

    val
end
