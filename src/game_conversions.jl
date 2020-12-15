"""
Facet (`Vector{Int64}`) -> `BellGame`

    convert(
        ::Type{BellGame},
        facet::Vector{Int64},
        scenario::Union{BlackBox, LocalSignaling};
        rep = "normalized"::String
    )
"""
function convert( ::Type{BellGame},
    facet::Vector{Int64},
    scenario::Union{BlackBox,LocalSignaling};
    rep = "normalized"::String
)
    if !(rep in ("normalized", "generalized"))
        throw(DomainError(rep, "Argument `rep` must be either 'normalized' or 'generalized'"))
    end

    game_dims = strategy_dims(scenario)

    div_facet = facet .÷ gcd(facet...)

    bound = div_facet[end]
    game_matrix = (rep == "normalized") ? cat(
        reshape(div_facet[1:end-1], (game_dims[1]-1, game_dims[2])),
        zeros(Int64, (1,game_dims[2])),
        dims=1
    ) : reshape(div_facet[1:end-1], game_dims)

    BellGame(_reverse_game_normalization(game_matrix, bound)...)
end

"""
Facet (`Vector{Int64}`) -> `BellGame`

    convert(
        ::Type{BellGame},
        facet::Vector{Int64},
        scenario::BipartiteNoSignaling;
        rep = "no-signaling"::String
    )

Transforms LocalPolytope facets into `BellGame`  types.
"""
function convert(::Type{BellGame},
    facet::Vector{Int64},
    scenario::BipartiteNoSignaling;
    rep = "no-signaling"::String
)
    if !(rep in ["no-signaling","normalized","generalized"])
        throw(DomainError(rep, "input `rep` must be in [\"no-signaling\",\"normalized\",\"generalized\"]"))
    end

    game_dims = strategy_dims(scenario)
    game = (rep == "generalized") ? reshape(facet[1:end-1], game_dims) : zeros(Int64, game_dims)
    bound = facet[end]

    if rep == "no-signaling"
        α_dim = (scenario.A-1)*scenario.X
        β_dim = (scenario.B-1)*scenario.Y

        α_game = reshape(facet[1:α_dim], (scenario.A-1, scenario.X))
        β_game = reshape(facet[α_dim+1:α_dim+β_dim], (scenario.B-1, scenario.Y))
        αβ_game = reshape(facet[α_dim+β_dim+1:end-1], ((scenario.A-1)*(scenario.B-1), scenario.X*scenario.Y))

        for i in 0:scenario.A-2
            game[i*scenario.B+1:i*scenario.B + scenario.B-1,:] = αβ_game[i*(scenario.B-1)+1:(i+1)*(scenario.B-1),:]

            b_vec = zeros(Int64, scenario.Y)
            b_vec[1] = 1

            game[i*scenario.B+1:(i+1)*scenario.B,:] += ones(Int64,scenario.B)*kron(α_game[i+1,:],b_vec)'
        end

        for i in 1:scenario.B-1
            a_vec = zeros(Int64, scenario.X)
            a_vec[1] = 1

            game_row_ids = i:scenario.B:scenario.A*scenario.B-1
            αβ_game_row_ids = i:scenario.B-1:(scenario.A-1)*(scenario.B-1)

            game[game_row_ids,:] += ones(Int64,scenario.A)*kron(a_vec,β_game[i,:])'
        end
    elseif rep == "normalized"
        game[1:game_dims[1]-1,:] = reshape(facet[1:end-1], (game_dims[1]-1, game_dims[2]))
    end

    BellGame(_reverse_game_normalization(game,bound)...)
end

"""
`XPORTA.IEQ` to `BellGame`'s

    convert(
        ::Type{Vector{BellGame}},
        ieq::IEQ,
        scenario::Union{BlackBox, LocalSignaling};
        rep = "normalized" :: String
    )
"""
function convert(::Type{Vector{BellGame}},
    ieq::IEQ,
    scenario::Union{BlackBox,LocalSignaling};
    rep = "normalized"::String
)
    inequalities = convert.(Int64, ieq.inequalities)
    map( row_id -> convert(BellGame, inequalities[row_id,:], scenario, rep=rep), 1:size(inequalities,1))
end

"""
`BellGame` -> Facet (`Vector{Int64}`)

    convert(::Type{Vector{Int64}}, BG::BellGame; rep = "normalized" :: String)
"""
function convert(::Type{Vector{Int64}}, BG::BellGame; rep = "normalized"::String)
    if !(rep in ("normalized", "generalized"))
        throw(DomainError(rep, "Argument `rep` must be either 'normalized' or 'generalized'"))
    end

    bound = BG.β
    game_matrix = BG.game[:,:]

    game_dims = size(game_matrix)

    if rep == "normalized"
        for col_id in 1:game_dims[2]
            col = game_matrix[:,col_id]
            if col[end] !== 0
                game_matrix[:,col_id] .-= col[end]
                bound -= col[end]
            end
        end

        game_matrix = game_matrix[1:(end-1),:]
    end

    cat(game_matrix[:], bound, dims=1)
end

"""
BellGame -> Vector{Int64}

    convert(::Type{Vector{Int64}},
        BG::BellGame,
        scenario::BipartiteNoSignaling;
        rep = "no-signaling" :: String
    )

Transforms a `BellGame` for a `BipartiteNoSignaling` scenario into a facet vector.
"""
function convert(::Type{Vector{Int64}},
    BG::BellGame,
    scenario::BipartiteNoSignaling;
    rep = "no-signaling" :: String
)
    if !(rep in ["no-signaling","normalized","generalized"])
        throw(DomainError(rep, "input `rep` must be in [\"no-signaling\",\"normalized\",\"generalized\"]"))
    end

    game_dims = size(BG)
    v_dim = LocalPolytope.vertex_dims(scenario, rep)

    facet = (rep == "generalized") ? vcat(BG[:], BG.β) : zeros(Int64, v_dim+1)

    if rep == "normalized"
        (game_matrix, bound) = _apply_game_normalization(BG[:,:], BG.β)

        facet = vcat(game_matrix[1:game_dims[1]-1,:][:], bound)
    elseif rep == "no-signaling"
        (game_matrix, bound) = _apply_game_normalization!(BG[:,:], BG.β)

        # construct G(a|x) and G(b|y)
        # in each column, subtract off from each Alice/Bob column the values excluded from the no-signaling
        α_game = zeros(Int64, (scenario.A-1, scenario.X))
        β_game = zeros(Int64, (scenario.B-1, scenario.Y))

        for a in 1:scenario.A-1
            target_row = a * scenario.B
            subtract_vals = game_matrix[target_row,:]

            a_dims = (a-1)*scenario.B +1: a * scenario.B

            game_matrix[a_dims,:] -= ones(Int64, scenario.B) * subtract_vals'

            α_game_rows = map(x -> begin
                x_dims = (x-1)*scenario.Y+1:x*scenario.Y

                sum(subtract_vals[x_dims])
            end, 1:scenario.X)

            α_game[a,:] = α_game_rows
        end

        for b in 1:scenario.B-1
            target_row = (scenario.A-1) * (scenario.B) + b
            subtract_vals = game_matrix[target_row,:]

            b_dims = b:scenario.B: scenario.A * scenario.B -1

            game_matrix[b_dims,:] -= ones(Int64, scenario.A) * subtract_vals'

            β_game_rows = map(y -> begin
                y_dims = y:scenario.Y: scenario.X*y

                sum(subtract_vals[y_dims])
            end, 1:scenario.Y)

            β_game[b,:] = β_game_rows
        end

        # All remaining terms are in the no-sig subspace and are taken as is
        αβ_game = zeros(Int64, ((scenario.A-1)*(scenario.B-1), scenario.X*scenario.Y))
        for a in 1:scenario.A-1
            αβ_game[(a-1)*(scenario.B-1) + 1:a*(scenario.B-1),:] = game_matrix[(a-1)*scenario.B+1:a*scenario.B-1,:]
        end

        facet = vcat(α_game[:], β_game[:], αβ_game[:], bound)
    end

    facet
end

"""
`BellGame`'s to `XPORTA.IEQ`

    convert(::Type{IEQ}, bell_games::Vector{BellGame}; rep = "normalized" :: String)
"""
function convert(::Type{IEQ}, bell_games::Vector{BellGame}; rep = "normalized"::String)
    ieq_vectors = map( bg -> convert(Vector{Int64}, bg, rep=rep), bell_games )
    IEQ(inequalities = hcat(ieq_vectors...)'[:,:])
end

"""

Applies the normalization constraint to remove all negative values in the provided
`game_matrix`. Returns a tuple `(new_game_matrix, new_bound)`
"""
function _reverse_game_normalization(game_matrix::Matrix{Int64}, bound::Int64) :: Tuple{Matrix{Int64}, Int64}
    new_bound = bound
    new_game_matrix = game_matrix

    for col_id in 1:size(game_matrix,2)
        col = game_matrix[:,col_id]
        col_min =  min(col...)

        if col_min != 0
            new_game_matrix[:,col_id] .-= col_min
            new_bound -= col_min
        end
    end

    (new_game_matrix, new_bound)
end

function _apply_game_normalization!(game_matrix::Matrix{Int64}, bound::Int64) :: Tuple{Matrix{Int64}, Int64}

    for col_id in 1:size(game_matrix,2)
        col = game_matrix[:,col_id]
        if col[end] !== 0
            game_matrix[:,col_id] .-= col[end]
            bound -= col[end]
        end
    end

    (game_matrix, bound)
end
