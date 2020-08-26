export  dimension

"""
    dimension( vertices :: Vector{Vector{Int64}} ) :: Int64

For the provided vertices (points), finds the dimension of the affine space spanned
by the vertices. This method computes the dimension of a polytope, facet, or
collection of points.

This also accepts matrices and `DeterministicStrategy` types as arguments:

    dimension( vertices :: Vector{Matrix{Int64}} ) :: Int64

    dimension( vertices :: Vector{DeterministicStrategy} ) :: Int64
"""
function dimension(vertices :: Vector{Vector{Int64}}) :: Int64
    vectors = map( v -> v - vertices[1], vertices[2:end])
    rank(hcat(vectors...))
end

function dimension(vertices :: Vector{Matrix{Int64}}) :: Int64
    dimension(map( v -> v[:], vertices))
end

function dimension(vertices :: Vector{DeterministicStrategy}) :: Int64
    dimension(map(v -> v[:], vertices))
end
