"""
"""
module Symmetry

using ..BellScenario: Scenario, PrepareAndMeasure, BellGame

using ..QBase: QMath

abstract type Group end


"""

Finds the canonical facet for the
"""
function generator_facet(BG :: BellGame, PM :: PrepareAndMeasure) :: BellGame

    game_copy = copy(BG)

    unique_vals = sort(unique(game_copy))
    index_maps = filter(tuple -> tuple[1] != tuple[2], map(i -> (unique_vals[i+1], i), 0:(length(unique_vals)-1)))

    # reducing values to play nice with lexicographic scoring
    for id_map in index_maps
        game_copy[ game_copy .== id_map[1] ] .= id_map[2]
    end

    target_row = 1
    col_perms = [collect(1:PM.X)]
    perm_game = _perm_increase_lexico_score(game_copy, target_row, col_perms)

    # correcting reduced values
    for id_map in reverse(index_maps)
        perm_game[ perm_game .== id_map[2] ] .= id_map[1]
    end

    BellGame(perm_game, BG.β)
end

function _perm_increase_lexico_score(
    game::Matrix{Int64},
    target_row::Int64,
    allowed_col_perms::Array{Array{Int64,1},1}
) :: Matrix{Int64}
    (num_rows, num_cols) =  size(game)
    base = max(game...) + 1

    game_row_dict = Dict()
    for row_id in target_row:num_rows
        row = game[row_id,:]
        if haskey(game_row_dict, row)
            push!(game_row_dict[row]["ids"], row_id)
        else
            game_row_dict[row] = Dict(
                "row_tuple" => map(col_id -> (game[row_id,col_id], col_id), 1:num_cols),
                "ids"  => [row_id]
            )
        end
    end

    # tag each row and row element with an id to track permutations.
    game_row_tuples = map(dict -> (dict["row_tuple"], dict["ids"]) , values(game_row_dict))

    # sort each row by reverse order to maximize lexicographic contribution
    sorted_game_row_tuples = map(game_row_tuple ->
            (_sort_perm_groups(game_row_tuple[1], allowed_col_perms), game_row_tuple[2]),
        game_row_tuples)

    # find the row which maximizes the lexicographic score and has the most duplicate rows
    max_game_row = sort(
            sort(sorted_game_row_tuples, by=tuple->length(tuple[2]), rev=true),
            rev = true,
            by = row_tuple -> QMath.base_n_val(map(tuple -> tuple[1], row_tuple[1]), base)
        )[1]

    # getting the row id to use for the target row
    max_game_row_ids = max_game_row[2]

    # getting the column permutations which maximize target row
    permute_ids = map(tuple -> tuple[2], max_game_row[1])

    # permuting columns and shifting maximal row to first row
    game[target_row:end,:] = game[[max_game_row_ids..., filter(x -> !in(x,max_game_row_ids), target_row:num_rows)...], permute_ids]

    # reducing the number of allowed column permutations
    col_perm_groups = Array{Array{Int64,1}}(undef,0)
    for perm_group in allowed_col_perms
        if length(perm_group) == 1
            # if we couldn't permute before we still can't
            push!(col_perm_groups, perm_group)
        elseif length(unique(game[target_row, perm_group])) == 1
            # if all elements are the same, we continue to  consider the perm group
            push!(col_perm_groups, perm_group)
        else
            for n in (base-1):-1:0
                sub_group_ids = findall( el -> el == n, game[target_row, perm_group])

                if length(sub_group_ids) >=  1
                    push!(col_perm_groups, perm_group[sub_group_ids])
                end
            end
        end
    end

    # recursive continuation condition, increase target riw and update permutations
    new_target_row = target_row + length(max_game_row_ids)
    if new_target_row <= num_rows
        game = _perm_increase_lexico_score(game, new_target_row, col_perm_groups)
    end

    game
end

function _sort_perm_groups(row_tuples::Array{Tuple{Int64,Int64},1}, perm_groups::Array{Array{Int64,1},1})
    for perm_ids in perm_groups
        if length(perm_ids) > 1
            row_tuples[perm_ids] = sort(row_tuples[perm_ids], rev=true, by= tuple->tuple[1])
        end
    end

    row_tuples
end

function groups(PM :: PrepareAndMeasure)
    input_relabels = QMath.permutation_maps(PM.X)
    output_relabels = QMath.permutation_maps(PM.B)

    Dict(
        "input" => input_relabels,
        "output" => output_relabels
    )
end

function lexico_score(BG :: BellGame) :: Vector{Int64}

    game_copy = copy(BG)
    unique_vals = sort(unique(game_copy))
    index_maps = filter(tuple -> tuple[1] != tuple[2], map(i -> (unique_vals[i+1], i), 0:(length(unique_vals)-1)))

    for map in index_maps
        game_copy[ game_copy .== map[1] ] .= map[2]
    end

    base = length(unique_vals)

    map(row -> QMath.base_n_val(row, base), eachrow(game_copy))
end


end
