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

    scalar = gcd(a_rot..., β_rot)

    HalfSpace(Int.(a_rot/scalar), Int(β_rot/scalar))
end

function adjacent_facets(F::HalfSpace, vertices)
    # vertices of facet F
    F_vertices = filter(v -> F.a' * v == F.β, vertices)

    # find the facets of Facet F, these subfacets are labeled G
    G_ieq = traf(POI(vertices = vcat(map(v -> v', F_vertices)...)))
    G_ineqs = convert.(Int, G_ieq.inequalities)

    G_facets = [ HalfSpace(G_ineqs[i,1:(end-1)], G_ineqs[i,end]) for i in 1:size(G_ineqs,1) ]

    # polytope vertices not in F index 1 is farthest from F
    xbar_vertices = sort(filter( v -> F.a' * v != F.β, vertices), by = v -> F.a' * v)

    adjacent_facets = []

    for G in G_facets
        for xbar_id in 1:length(xbar_vertices)
            # find v in P such that v is is farthest from the facet (first in above list)
            # this vertex will be called x_bar
            xbar = xbar_vertices[xbar_id]

            # rotate G about x_bar to get G'
            G_rot = rotate_facet(F, G, xbar)

            # sort xbars by score against facet, greatest to least
            sorted_xbars = sort(xbar_vertices, by = v -> G_rot.a' * v, rev=true)

            println(G_rot.a' * sorted_xbars[1], " is this the bound ", G_rot.β)

            # if max xbar is on the facet, then G_rot is a polytope facet adjacent to F
            if G_rot.a' * sorted_xbars[1] == G_rot.β
                push!(adjacent_facets, G_rot)
                break
            end
        end
    end


    return Dict(
        "facet" => F,
        "facet_vertices" => F_vertices,
        "xbar_vertices" => xbar_vertices,
        "G_facets" => G_facets,
        "adjacent_facets" => adjacent_facets,
    )
end

function adjacency_decomposition(facet, vertices)

end
