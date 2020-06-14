using XPORTA, Polyhedra


"""
    rotate_facet(F::HalfSpace, G::HalfSpace xbar::Vector{Int})

Performs a rotation of facet `F` relative to non-included vertex `xbar`. Returns
the rotated `HalfSpace`. `F` is a polytope facet, `G` is a facet of `F` and `xbar`
is a vertex not contained on the boundary of `F`.
"""
function rotate_facet(F::HalfSpace, G::HalfSpace, xbar::Vector{Int}) :: HalfSpace
    a_rot = (F.β - F.a'*xbar)*G.a - (G.β - G.a'*xbar)*F.a
    β_rot = (F.β - F.a'*xbar)*G.β - (G.β - G.a'*xbar)*F.β

    HalfSpace(a_rot, β_rot)
end

function adjacent_facets(facet, vertices)

end

function adjacency_decomposition(facet, vertices)

end
