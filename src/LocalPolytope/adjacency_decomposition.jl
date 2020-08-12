using XPORTA

export rotate_facet, adjacent_facets, adjacency_decomposition

"""
    rotate_facet(
        F :: Vector{Int64},
        G :: Vector{Int64},
        xbar :: Vector{Int64}
    ) :: Vector{Int64}

Performs a rotation of facet `F` relative to non-included vertex `xbar` and returns
the rotated facet vector. `F` is a polytope facet, `G` is a subfacet of `F` and `xbar`
is a vertex not contained by `F`. By construction, the rotated facet contains
`xbar`.
"""
function rotate_facet(F::Vector{Int64}, G ::Vector{Int64}, xbar::Vector{Int64}) :: Vector{Int64}
    a_rot = (F[end] - F[1:(end-1)]' * xbar)*G[1:(end-1)] - (G[end] - G[1:(end-1)]' * xbar)*F[1:(end-1)]
    β_rot = (F[end] - F[1:(end-1)]' * xbar)*G[end] - (G[end] - G[1:(end-1)]' * xbar)*F[end]

    cat(a_rot, β_rot, dims=1)
end

"""
    adjacent_facets(
        vertices :: Vector{Vector{Int64}},
        F :: Vector{int64},
    ) :: Vector{Vector{Int64}}

For the polytope represented by `vertices`, returns the set of facets
adjacent to `F`.

Facet vector `F` and the return facet vectors are assumed to ba in the normalized
subspace.
"""
function adjacent_facets(vertices::Vector{Vector{Int64}}, F::Vector{Int64}) :: Vector{Vector{Int64}}
    # vertices of facet F
    F_vertices = filter(v -> F[1:(end-1)]' * v == F[end], vertices)

    # find the subfacets of facet F, these subfacets are labeled G
    G_ieq = traf(POI(vertices = hcat(F_vertices...)'[:,:]))
    G_ineqs = convert.(Int64, G_ieq.inequalities)

    # polytope vertices not in F index 1 is farthest from F
    xbar_vertices = sort(
        filter(v -> F[1:(end-1)]' * v != F[end], vertices),
        by = v -> F[1:(end-1)]' * v, rev=true
    )

    adjacent_facets = Array{Vector{Int64},1}(undef,0)
    for G_row_id in 1:size(G_ineqs,1)
        G = G_ineqs[G_row_id,:]
        for xbar_id in 1:length(xbar_vertices)
            xbar = xbar_vertices[xbar_id]

            # rotate G about x_bar to get G_rot
            G_rot = rotate_facet(F, G, xbar)

            # If no polytope vertices violate G_rot facet then it is a facet of polytope P
            if (findfirst(v -> G_rot[1:(end-1)]' * v > G_rot[end], xbar_vertices) === nothing)
                push!(adjacent_facets, G_rot)

                # breaking misses adjacent facets, but I have not found any cases
                # where a canonical facet class has been missed
                break # TODO: does breaking catch all adjacent facet classes?
            end
        end
    end

    return adjacent_facets
end

"""
    adjacenecy_decomposition(
        vertices :: Vector{Vector{Int64}},
        BG_seed :: BellGame,
        PM :: PrepareAndMeasure;
        kwargs
    )

Given a polytpe represented by `vertices`, returns the complete set of canonical
facets for prepare and measure scenario `PM`. The adjacency_decomposition algorithm
requires a seeded vertex which is supplied with the `BG_seed` argument. Facets
are returned in the lexicographic normal form.

Returns a dictionary where the keys are canonical `BellGames` and the value is a
dictionary with keys

* "considered" => true, if the facet was considered in the algorithm.
* "skipped" => true, if the facet was skipped (not considered).
* "num_vertices" => number of vertices.
* "norm_facet" => facet vector representative (normalized rep) of canonical game

Keyword  arguments `kwargs`
* `skip_games ::  Vector{BellGame}` - Optional list of games to skip.
* `max_vertices :: Int64` - Defaults to 100, the maximum number of vertices to allow in target facets.
"""
function adjacency_decomposition(
    vertices :: Vector{Vector{Int64}},
    BG_seed :: BellGame,
    PM :: PrepareAndMeasure;
    skip_games = Array{BellGame}(undef,0) :: Vector{BellGame},
    max_vertices = 100 :: Int64
)
    # canonicalize facet
    canonical_BG_seed = LocalPolytope.generator_facet(BG_seed, PM)
    canonical_skip_BGs = map( BG -> LocalPolytope.generator_facet(BG, PM), skip_games)

    norm_facet_seed = convert(Vector{Int64}, canonical_BG_seed, rep = "normalized")
    num_BG_seed_vertices = length(
        filter(v -> norm_facet_seed[1:(end-1)]' * v == norm_facet_seed[end], vertices)
    )

    # holds considered and unconsidered facets
    facet_dict = Dict{BellGame, Dict}(
        canonical_BG_seed => Dict(
            "considered" => false,
            "skipped" => false,
            "norm_facet"  => norm_facet_seed,
            "num_vertices" => num_BG_seed_vertices
        )
    )

    # add skipped facets as considered
    for skip_BG in canonical_skip_BGs
        facet_dict[skip_BG] = Dict("considered" => true, "skipped" => true)
    end

    # loop until all facet are considered
    while !all(d -> d["considered"], collect(values(facet_dict)))
        # select the first uncosidered facet with fewest vertices
        target_BG = sort(
                filter(d -> facet_dict[d]["considered"] == false, collect(keys(facet_dict))),
                by=key -> facet_dict[key]["num_vertices"]
            )[1]

        # get target facet vector in normalized representation
        norm_facet = facet_dict[target_BG]["norm_facet"]

        # compute adjacent facets
        adj_facets = adjacent_facets(vertices, norm_facet)

        for adj_facet in adj_facets
            adj_game = convert(BellGame, adj_facet, PM, rep="normalized")
            canonical_game = LocalPolytope.generator_facet(adj_game, PM)

            # add new facet to  dictionary
            if !haskey(facet_dict, canonical_game)
                num_vertices = length(filter(v -> adj_facet[1:(end-1)]' * v == adj_facet[end], vertices))

                if num_vertices <= max_vertices
                    facet_dict[canonical_game] = Dict(
                        "considered" => false,
                        "skipped" => false,
                        "norm_facet" => adj_facet,
                        "num_vertices" => num_vertices
                    )
                else
                    facet_dict[new_game] = Dict(
                        "considered" => true,
                        "skipped" => true,
                        "norm_facet" => adj_facet,
                        "num_vertices" => num_vertices
                    )
                end
            end
        end

        # set current bell game "considered" to true
        facet_dict[target_BG]["considered"] = true
    end

    facet_dict
end
