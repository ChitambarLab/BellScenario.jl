using XPORTA: POI, IEQ, traf, make_porta_tmp, rm_porta_tmp
using JSON, Dates

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
        F :: Vector{int64};
        dir = "./" :: String,
        cleanup = true ::Bool
    ) :: Vector{Vector{Int64}}

For the polytope represented by `vertices`, returns the set of facets
adjacent to `F`.

Facet vector `F` and the return facet vectors are assumed to ba in the normalized
subspace.

The  `dir` argument specifies where to where to write files and directories from
XPORTA.jl. If `cleanup` is `true`, then a  porta_tmp directory is created as a
subdirectory of `dir`.

If `cleanup` is `false`, the created  `porta_tmp` directory is not removed.
"""
function adjacent_facets(
    vertices::Vector{Vector{Int64}},
    F::Vector{Int64};
    dir = "./" :: String,
    cleanup=true :: Bool
) :: Vector{Vector{Int64}}
    # vertices of facet F
    F_vertices = filter(v -> F[1:(end-1)]' * v == F[end], vertices)

    println("befor traf")


    # find the subfacets of facet F, these subfacets are labeled G
    G_ieq = traf(POI(vertices = hcat(F_vertices...)'[:,:]), dir=dir, cleanup=cleanup)
    G_ineqs = convert.(Int64, G_ieq.inequalities)

    println("after traf")

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
        scenario :: LocalSignaling;
        kwargs...
    ) :: Dict

Given a polytpe represented by `vertices`, returns the complete set of canonical
facets for prepare and measure scenario `scenario`. The adjacency_decomposition algorithm
requires a seeded vertex which is supplied with the `BG_seed` argument. Facets
are returned in the lexicographic normal form.

### Returned Dictionary Format

Returns a dictionary where the keys are canonical `BellGames` and the value is a
dictionary with keys

* "considered" => true, if the facet was considered in the algorithm.
* "skipped" => true, if the facet was skipped (not considered).
* "num_vertices" => number of vertices.
* "norm_facet" => facet vector representative (normalized rep) of canonical game

### Keyword  Arguments: `kwargs...`
* `skip_games = [] ::  Vector{BellGame}` - List of games to skip.
* `max_vertices = 100 :: Int64` - The maximum number of vertices to allow in target facets.
* `dir` = "./" :: String` - Directory in which to create `porta_tmp/` and `.json` files.
* `log = false :: Bool` - If true, the facet dictionary is  written to `.json` each iteration.
* `log_filename = "adjacency_decomposition_$(Dates.now).json" :: String`
"""
function adjacency_decomposition(
    vertices :: Vector{Vector{Int64}},
    BG_seed :: BellGame,
    scenario :: LocalSignaling;
    skip_games = Array{BellGame}(undef,0) :: Vector{BellGame},
    max_vertices = 100 :: Int64,
    dir = "./"  :: String,
    log = false :: Bool,
    log_filename = "adjacency_decomposition_$(Dates.now()).json" :: String,
) :: Dict
    # canonicalize facet
    canonical_BG_seed = LocalPolytope.generator_facet(BG_seed, scenario)
    canonical_skip_BGs = map( BG -> LocalPolytope.generator_facet(BG, scenario), skip_games)

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
            "generator_facet" => convert(Vector{Int64}, canonical_BG_seed),
            "num_vertices" => num_BG_seed_vertices
        )
    )

    # add skipped facets as considered
    for skip_BG in canonical_skip_BGs
        facet_dict[skip_BG] = Dict(
            "considered" => true,
            "norm_facet" => convert(Vector{Int64}, skip_BG, rep = "normalized"),
            "generator_facet" => convert(Vector{Int64}, skip_BG),
            "skipped" => true
        )
    end

    # create porta_tmp directory
    porta_tmp_dir = make_porta_tmp(dir)

    # loop until all facet are considered
    while !all(d -> d["considered"], collect(values(facet_dict)))
        println("starting loop")

        # select the first uncosidered facet with fewest vertices
        target_BG = sort(
                filter(d -> facet_dict[d]["considered"] == false, collect(keys(facet_dict))),
                by=key -> facet_dict[key]["num_vertices"]
            )[1]

        # get target facet vector in normalized representation
        norm_facet = facet_dict[target_BG]["norm_facet"]

        println("about to enter try block")

        # compute adjacent facets
        adj_facets = try
            println("inside try")
            adjacent_facets(vertices, norm_facet, dir=porta_tmp_dir, cleanup=false)
        # if an unexpected error occurs with XPORTA, mark facet as such and move on.
        catch error
            facet_dict[target_BG]["considered"] = true

            push!(facet_dict[target_BG],"error" => true)
            push!(facet_dict[target_BG],"error_msg" => sprint(showerror, error, backtrace()))

            continue
        end

        for adj_facet in adj_facets
            adj_game = convert(BellGame, adj_facet, scenario, rep="normalized")
            canonical_game = LocalPolytope.generator_facet(adj_game, scenario)

            # add new facet to  dictionary
            if !haskey(facet_dict, canonical_game)
                num_vertices = length(filter(v -> adj_facet[1:(end-1)]' * v == adj_facet[end], vertices))

                if num_vertices <= max_vertices
                    facet_dict[canonical_game] = Dict(
                        "considered" => false,
                        "skipped" => false,
                        "norm_facet" => adj_facet,
                        "generator_facet" => convert(Vector{Int64}, canonical_game),
                        "num_vertices" => num_vertices
                    )
                else
                    facet_dict[canonical_game] = Dict(
                        "considered" => true,
                        "skipped" => true,
                        "norm_facet" => adj_facet,
                        "generator_facet" => convert(Vector{Int64}, canonical_game),
                        "num_vertices" => num_vertices
                    )
                end
            end
        end

        # set current bell game "considered" to true
        facet_dict[target_BG]["considered"] = true

        if log
            open(dir*log_filename, "w") do io
                JSON.print(io, facet_dict)
            end
        end
    end

    # cleanup porta_tmp directory after completion
    rm_porta_tmp(dir)

    facet_dict
end
