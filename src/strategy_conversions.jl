"""
Vertex (`Vector{Int64}`) -> `DeterministicStrategy`

    convert(
        ::Type{DeterministicStrategy},
        vertex  :: Vector{Int64},
        scenario :: Scenario;
        rep = "normalized" :: String
    )
"""
function convert(::Type{DeterministicStrategy}, vertex::Vector{Int64}, scenario::Union{LocalSignaling,BlackBox}; rep="normalized" :: String)
    if !(rep in ("normalized", "generalized"))
        throw(DomainError(rep, "Argument `rep` must be either 'normalized' or 'generalized'"))
    end

    s_dims = strategy_dims(scenario)

    strategy_matrix = (rep == "normalized") ? cat(
        reshape(vertex, (s_dims[1]-1, s_dims[2])),
        zeros(Int64, (1,s_dims[2])),
        dims = 1
    ) : reshape(vertex, s_dims)

    # for vertices in the normalized representation, if a column does not contain
    # a one, add it to the last row.
    if rep == "normalized"
        for col_id in 1:s_dims[2]
            if sum(strategy_matrix[:,col_id]) == 0
                strategy_matrix[s_dims[1], col_id] = 1
            end
        end
    end

    DeterministicStrategy(strategy_matrix, scenario)
end

"""
`DeterministicStrategy` -> Vertex (`Vector{Int64}`)

    convert(
        ::Type{Vector{Int64}},
        strategy :: DeterministicStrategy;
        rep = "normalized" :: String
    )
"""
function convert(::Type{Vector{Int64}}, strategy::DeterministicStrategy; rep = "normalized"::String)
    if !(rep in ("normalized", "generalized"))
        throw(DomainError(rep, "Argument `rep` must be either 'normalized' or 'generalized'"))
    end

    max_row = (rep == "normalized") ? size(strategy,1) - 1 : size(strategy,1)

    strategy[1:max_row,:][:]
end

"""
    convert(
        S :: Type{<:AbstractStrategy}, vertex::Vector{<:Real}, scenario::BipartiteNonSignaling;
        rep="non-signaling" :: String
    )

Transforms a behavior vector or vertex in to either a [`Strategy`](@ref) or [`DeterministicStrategy`](@ref).
If converting into a `DeterministicStrategy`, the vertex must contain `Int64` values.
Valid representations are "non-signaling", "normalized", and "generalized".
"""
function convert(S::Type{<:AbstractStrategy}, vertex::Vector{<:Real}, scenario::BipartiteNonSignaling;
    rep="non-signaling" :: String
)
    if !(rep in ["non-signaling","normalized","generalized"])
        throw(DomainError(rep, "input `rep` must be in [\"non-signaling\",\"normalized\",\"generalized\"]"))
    end

    num_type = (S == Strategy) ? Float64 : Int64

    if S == DeterministicStrategy && eltype(vertex) !== num_type
        throw(DomainError((num_type), "for type $S, vertex must be `Vector{$num_type}`."))
    end

    gen_strat_dims = strategy_dims(scenario)
    gen_strat = (rep == "generalized") ? reshape(vertex, gen_strat_dims) : zeros(num_type, gen_strat_dims)

    if rep == "non-signaling"
        # code to get non-signaling to normalized
        α_dim = (scenario.A-1)*scenario.X
        β_dim = (scenario.B-1)*scenario.Y

        α_strat = reshape(vertex[1:α_dim], (scenario.A-1, scenario.X))
        β_strat = reshape(vertex[α_dim+1:α_dim+β_dim], (scenario.B-1, scenario.Y))
        αβ_strat = reshape(vertex[α_dim+β_dim+1:end], ((scenario.A-1)*(scenario.B-1), scenario.X*scenario.Y))

        # reversing non-signaling constraints for Alice
        for i in 0:scenario.A-2
            gen_strat[i*scenario.B+1:i*scenario.B + scenario.B-1,:] = αβ_strat[i*(scenario.B-1)+1:(i+1)*(scenario.B-1),:]

            p_ax_row = kron(α_strat[i+1,:],ones(Int64, scenario.Y))'
            gen_strat[i*scenario.B+scenario.B,:] = p_ax_row - ones(Int64, (1,scenario.B-1)) * αβ_strat[i*(scenario.B-1)+1:(i+1)*(scenario.B-1),:]
        end

        for i in 1:scenario.B-1
            αβ_row_ids = i:scenario.B-1:(scenario.A-1)*(scenario.B-1)

            p_by_row = kron(ones(Int64,scenario.X),β_strat[i,:])'
            gen_strat[(scenario.A-1)*scenario.B + i,:] = p_by_row - ones(Int64, (1,scenario.A-1)) * αβ_strat[αβ_row_ids,:]
        end

    elseif rep == "normalized"
        gen_strat[1:gen_strat_dims[1]-1,:] = reshape(vertex, (gen_strat_dims[1]-1, gen_strat_dims[2]))
    end

    # strategy is in "normalized" rep transform it to "generalized" rep
    if rep in ("normalized", "non-signaling")
        gen_strat[end,:] = 1 .- ones(Int64, (1,gen_strat_dims[1]-1)) * gen_strat[1:end-1,:  ]
    end

    (S == Strategy) ? Strategy(gen_strat) : DeterministicStrategy(gen_strat)
end
