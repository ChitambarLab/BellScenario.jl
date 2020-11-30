export facets

"""
    facets(
        vertices :: Vector{Vector{Int64}};
        dir = "./" :: String,
        cleanup=true :: Bool
    ) :: Dict{String, Vector{Vector{Int64}}}

Computes the facet inequalities and equalities which bound the convex polyhedron
defined by the set of `vertices`.
If the optimal representation is used for the vertices, then no equalitities will
be present.
For example, if the "generalized" vertex representation is used the normalization
constraints will be included in the resulting equalities.
The output of this method is structured:

```
Dict(
    "facets" => inequalities, # :: Vector{Vector{Int64}}
    "equalitites" => equalities, # :: Vector{Vector{Int64}}
)
```

Facet inequalities and equalities are represented as single vector `f` where
`f[1:(end-1)]` contains the coefficients of the linear inquality and `f[end]`
is the bound.
Facet inequalities and equalities act upon a vertex `v` where inequalities are
arranged such that there is an implicit `f[1:(end-1)]' * v ≤ f[end]`
and equalities are arranged such that `f[1:(end-1)]' * v == f[end]`

!!! note "Supporting Software"
    The vertex -> facet transformation is performed using the
    [`traf`](https://juliapolyhedra.github.io/XPORTA.jl/dev/exports/#XPORTA.traf)
    method of [XPORTA.jl](https://github.com/JuliaPolyhedra/XPORTA.jl).
    Please refer to the source code for more details.

!!! warning "Performance Limitations"
    Vertex -> facet transformations are notoriously difficult computational problems.
    The performance limits of this method will be met far before the limits the
    vertex enumeration methods.
"""
function facets(
    vertices :: Vector{Vector{Int64}};
    dir = "./" :: String,
    cleanup=true :: Bool
) :: Dict{String, Vector{Vector{Int64}}}
    ieq = traf(POI(vertices = hcat(vertices...)'[:,:]), dir=dir, cleanup=cleanup)

    inequalities = convert.(Int64, ieq.inequalities)
    equalities = convert.(Int64, ieq.equalities)

    Dict(
        "facets" => map(row_id -> inequalities[row_id,:], 1:size(inequalities,1)),
        "equalities" => map(row_id -> equalities[row_id,:], 1:size(equalities,1))
    )
end
