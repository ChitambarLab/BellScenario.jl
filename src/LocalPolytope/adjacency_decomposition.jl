using XPORTA

using Polyhedra: HalfSpace

using ..BellScenario: Scenario, PrepareAndMeasure, BellGame,  Behavior

"""
    rotate_facet(F::HalfSpace, G::HalfSpace, xbar::Vector{Int})

Performs a rotation of facet `F` relative to non-included vertex `xbar`. Returns
the rotated `HalfSpace`. `F` is a polytope facet, `G` is a facet of `F` and `xbar`
is a vertex not contained on the boundary of `F`.
"""
function rotate_facet(F::HalfSpace, G::HalfSpace, xbar::Vector{Int}) :: HalfSpace
    a_rot = (F.β - F.a'*xbar)*G.a - (G.β - G.a'*xbar)*F.a
    β_rot = (F.β - F.a'*xbar)*G.β - (G.β - G.a'*xbar)*F.β

    HalfSpace(Int.(a_rot), Int(β_rot))
end

"""
    adjacent_facets( vertices, F::HalfSpace, PM::PrepareAndMeasure )

For the polytope represented by `vertices`, returns the canonical set of facets
adjacent to `F` for the specified prepare and measure scenario `PM`.

`F` is assumed to ba a valid facet in the normalized subspace.
"""
function adjacent_facets(vertices, F::HalfSpace, PM::PrepareAndMeasure)
    # vertices of facet F
    F_vertices = filter(v -> (F.a' * v)[1] == F.β, vertices)

    # find the facets of Facet F, these subfacets are labeled G
    G_ieq = traf(POI(vertices = vcat(map(v -> v', F_vertices)...)))
    G_ineqs = convert.(Int, G_ieq.inequalities)

    G_facets = [ HalfSpace(G_ineqs[i,1:(end-1)], G_ineqs[i,end]) for i in 1:size(G_ineqs,1) ]

    # polytope vertices not in F index 1 is farthest from F
    xbar_vertices = sort(filter( v -> (F.a' * v)[1] != F.β, vertices), by = v -> (F.a' * v)[1], rev=true)

    adjacent_facets = Array{BellGame,1}(undef,0)

    for G in G_facets
        for xbar_id in 1:length(xbar_vertices)
            # find v in P such that v is is farthest from the facet (first in above list)
            # this vertex will be called x_bar
            xbar = convert.(Int,xbar_vertices[xbar_id][:])

            # rotate G about x_bar to get G'
            G_rot = rotate_facet(F, G, xbar)

            # If no polytope vertices violate G_rot facet then it is a facet of polytope P
            if (findfirst(v -> (G_rot.a'*v)[1] > G_rot.β, xbar_vertices) === nothing)

                scalar = gcd(G_rot.a..., G_rot.β)

                gen_facet = cat([-1*G_rot.β/scalar], G_rot.a/scalar, zeros(Int64, PM.X), dims=1)

                (bound, facet_matrix) = LocalPolytope.facet_to_matrix(PM.X, PM.B, gen_facet)

                canonical_facet = LocalPolytope.generator_facet(BellGame(convert.(Int64,facet_matrix), convert(Int64,bound)), PM)

                push!(adjacent_facets, canonical_facet)
                break
            end
        end
    end

    # TODO: proactively handle unique facets
    return unique(adjacent_facets)
end

"""
    adjacenecy_decomposition(vertices, BG_seed::BellGame, PM::PrepareAndMeasure; kwargs )

Given a polytpe represented by `vertices`, returns the complete set of canonical
facets for prepare and measure scenario `PM`. The adjacency_decomposition algorithm
requires a seeded vertex which is supplied with the `BG_seed` argument. Facets
are returned in the lexicographic normal form.

Keyword  arguments `kwargs`
* `skip_games ::  Vector{BellGame}` - Optional list of games to skip.
* `max_vertices :: Int64` - Defaults to 100, the maximum number of vertices to allow in target facets.
"""
function adjacency_decomposition(vertices, BG_seed::BellGame, PM::PrepareAndMeasure; skip_games=Array{BellGame}(undef,0)::Vector{BellGame}, max_vertices=100::Int64)
    # projects generalized facets to normalized facets
    facet_proj = Behavior.norm_to_gen_proj((PM.X,1),(1,PM.B))

    # canonicalize facet
    canonical_BG_seed = LocalPolytope.generator_facet(BG_seed, PM)
    canonical_skip_BGs = map( BG -> LocalPolytope.generator_facet(BG, PM), skip_games)

    gen_facet_seed = cat([-1*canonical_BG_seed.β],canonical_BG_seed'[:], dims=1)
    norm_facet_seed = gen_facet_seed'*facet_proj
    num_BG_seed_vertices = length(filter(v -> (norm_facet_seed[2:end]'*v)[1] == -1*norm_facet_seed[1], vertices))

    # holds considered and unconsidered facets
    facet_dict = Dict{BellGame, Dict}(canonical_BG_seed => Dict("considered" => false, "num_vertices" =>  num_BG_seed_vertices))

    # add skipped facets as considered
    for skip_BG in canonical_skip_BGs
        facet_dict[skip_BG] = Dict("considered" => true)
    end

    # loop until all facet are considered
    while !all(d -> d["considered"], collect(values(facet_dict)))
        # select the first uncosidered facet with fewest vertices
        target_BG = sort(
                filter(d -> facet_dict[d]["considered"] == false, collect(keys(facet_dict))),
                by=key -> facet_dict[key]["num_vertices"]
            )[1]

        # convert Bell Game to HalfSpace (generalized -> normalized rep)
        gen_facet = cat([-1*target_BG.β],target_BG'[:], dims=1)
        norm_facet = gen_facet'*facet_proj

        target_HS = HalfSpace(norm_facet[2:end], -1*norm_facet[1])

        # compute adjacent facets
        adj_games = adjacent_facets(vertices, target_HS, PM)

        # add new facets to dictionary
        new_games = filter(game -> !haskey(facet_dict, game), adj_games)
        for new_game in new_games
            new_gen_facet = cat([-1*new_game.β],new_game'[:], dims=1)
            new_norm_facet = new_gen_facet'*facet_proj
            num_vertices =  length(filter(v -> (new_norm_facet[2:end]'*v)[1] == -1*new_norm_facet[1], vertices))

            if num_vertices <= max_vertices
                facet_dict[new_game] = Dict("considered" => false,  "num_vertices" => num_vertices)
            else
                facet_dict[new_game] = Dict("considered" => true,  "num_vertices" => num_vertices)
            end
        end

        # set current bell game to true
        facet_dict[target_BG]["considered"] = true
    end

    collect(keys(facet_dict))
end
