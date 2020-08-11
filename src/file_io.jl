export pretty_print_txt

"""
    pretty_print_txt( bell_games :: Vector{BellGame}, filename :: String )

Prints a set of `BellGame`'s to `filename.txt` in a human-readable form.
"""
function pretty_print_txt(bell_games :: Vector{BellGame}, filename :: String)
    ext = occursin(r"\.txt$", filename) ? "" : ".txt"
    filepath = filename * ext
    open(filepath, "w") do io
        for id in 1:length(bell_games)
            game = bell_games[id]
            print(io, "(", id ,")  ", game.β, " ≥ ")
            print(io, "\n")

            for game_row in eachrow(game)
                println(io, "\t", game_row)
            end

            print(io, "\n")
        end
    end
end
