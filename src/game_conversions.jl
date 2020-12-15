"""
Facet (`Vector{Int64}`) -> `BellGame`

    convert(
        ::Type{BellGame},
        facet::Vector{Int64},
        scenario::Scenario;
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


    for col_id in 1:game_dims[2]
        col = game_matrix[:,col_id]
        col_min =  min(col...)

        if col_min != 0
            game_matrix[:,col_id] .-= col_min
            bound -= col_min
        end
    end

    BellGame(game_matrix, bound)
end

"""
`XPORTA.IEQ` to `BellGame`'s

    convert(
        ::Type{Vector{BellGame}},
        ieq::IEQ,
        scenario::Union{BlackBox,LocalSignaling};
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
`BellGame`'s to `XPORTA.IEQ`

    convert(::Type{IEQ}, bell_games::Vector{BellGame}; rep = "normalized" :: String)
"""
function convert(::Type{IEQ}, bell_games::Vector{BellGame}; rep = "normalized"::String)
    ieq_vectors = map( bg -> convert(Vector{Int64}, bg, rep=rep), bell_games )
    IEQ(inequalities = hcat(ieq_vectors...)'[:,:])
end
