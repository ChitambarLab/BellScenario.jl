using XPORTA, Polyhedra

using ..BellScenario: Scenario, PrepareAndMeasure, BellGame

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
    adjacent_facets(F::HalfSpace, vertices, PM::PrepareAndMeasure)

assumes vertices and facets are in the normalized subspace
"""
function adjacent_facets(F::HalfSpace, vertices, PM::PrepareAndMeasure)
    # vertices of facet F
    F_vertices = filter(v -> (F.a' * v)[1] == F.β, vertices)

    # find the facets of Facet F, these subfacets are labeled G
    G_ieq = traf(POI(vertices = vcat(map(v -> v', F_vertices)...)))
    G_ineqs = convert.(Int, G_ieq.inequalities)

    G_facets = [ HalfSpace(G_ineqs[i,1:(end-1)], G_ineqs[i,end]) for i in 1:size(G_ineqs,1) ]

    # polytope vertices not in F index 1 is farthest from F
    xbar_vertices = sort(filter( v -> (F.a' * v)[1] != F.β, vertices), by = v -> (F.a' * v)[1], rev=true)

    adjacent_facets = Array{BellGame,1}(undef,0)

    iteration = 0
    for G in G_facets
        iteration = iteration + 1
        println("subfacet iteration : ", iteration, "/",length(G_facets) )
        for xbar_id in 1:length(xbar_vertices)
            # find v in P such that v is is farthest from the facet (first in above list)
            # this vertex will be called x_bar
            xbar = convert.(Int,xbar_vertices[xbar_id][:])

            # rotate G about x_bar to get G'
            G_rot = rotate_facet(F, G, xbar)

            # if max xbar is on the facet, then G_rot is a polytope facet adjacent to F
            if (findfirst(v -> (G_rot.a'*v)[1] > G_rot.β, xbar_vertices) === nothing)
                # check that vertex lies  on facet
                max_id = findfirst(v -> (G_rot.a'*v)[1] == G_rot.β,xbar_vertices)
                if max_id !== nothing#
                    println("xbar_id : ",xbar_id)
                    println("amx_id : ", max_id)

                    scalar = gcd(G_rot.a..., G_rot.β)

                    gen_facet = cat([-1*G_rot.β/scalar], G_rot.a/scalar, zeros(Int64, PM.X), dims=1)

                    (bound, facet_matrix) = LocalPolytope.facet_to_matrix(PM.X, PM.B, gen_facet)

                    canonical_facet = Symmetry.generator_facet(BellGame(convert.(Int64,facet_matrix), convert(Int64,bound)), PM)

                    push!(adjacent_facets, canonical_facet)
                    break
                end
            end
        end
    end

    return adjacent_facets
end

function adjacency_decomposition(BG_seed::BellGame, vertices, PM::PrepareAndMeasure)

    # canonicalize facet
    canonical_BG_seed  = Symmetry.generator_facet(BG_seed, PM)

    # holds considered facets
    facet_dict = Dict{BellGame, Bool}(canonical_BG_seed => false)

    # projects generalized facets to normalized facets
    facet_proj = Behavior.norm_to_gen_proj((PM.X,1),(1,PM.B))

    positivity_facet = zeros(Int64, (PM.B, PM.X))
    positivity_facet[1:(end-1), 1] .= 1

    positivity_game = BellGame( positivity_facet , 1)

    # loop until all facet are considered
    iteration = 0
    brute_count = 0
    while !all(values(facet_dict))
        # select the first uncosidered facet
        target_BG = findfirst(isequal(false), facet_dict)

        println("TARGET BELL GAME : ", target_BG, ", val : ", target_BG.β)

        # convert Bell Game to HalfSpace (generalized -> normalized rep)
        gen_facet = cat([-1*target_BG.β],target_BG'[:], dims=1)
        norm_facet = gen_facet'*facet_proj

        target_HS = HalfSpace(norm_facet[2:end], -1*norm_facet[1])

        # compute adjacent facets
        adj_games = adjacent_facets(target_HS, vertices, PM)

        # add new facets to dictionary
        new_games = filter(game -> !haskey(facet_dict, game), adj_games)
        for new_game in new_games

            if new_game != positivity_game
                facet_dict[new_game] = false
            end
        end

        # set current bell game to true
        facet_dict[target_BG] = true

        iteration = iteration + 1
        println("iteration : ", iteration)
        println("num_facets complete : ", length(filter(isequal(true), collect(values(facet_dict)))))
        println("num facets to go: ", length(filter(isequal(false), collect(values(facet_dict)))))
    end

    collect(keys(facet_dict))

end
