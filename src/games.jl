using XPORTA: IEQ

export AbstractGame, Game, BellGame

"""
Games represent a linear inequality.
"""
abstract type AbstractGame{T} <: AbstractMatrix{T} end
Base.size(G::AbstractGame) = size(G.game)
Base.getindex(G::AbstractGame, I::Vararg{Int64,2}) = getindex(G.game, I...)
Base.setindex!(G::AbstractGame, v, I::Vararg{Int64,2}) = (G.game[I...] = v)

"""
    Game(game::Matrix{T}, β::Real) <: AbstractMatrix{T}

A `Game` is played by a `Matrix` and `β` is the bounding or winning score of the
game. A `Game` is matrix representation of a linear inequality. Each element is
a linear scale factor for an element of a strategy matrix.

Type parameter `T` is typically an `Int64` or `Float64`.
"""
struct Game{T} <: AbstractGame{T}
    game :: Matrix{T}
    β :: Real
    scenario :: Scenario
    Game(game::Matrix{T}, β::Real) where T = new{T}(game, β)
end

"""
    BellGame(game::Matrix{Int64}, β::Int64, scenario::Scenario)

A tight Bell inequality for the correlation polytope for `scenario`.
"""
struct BellGame <: AbstractGame{Int64}
    game :: Matrix{Int64}
    β :: Int64
    scenario :: Scenario
    BellGame(game::Matrix{Int64}, β::Int64) = new(game, β)
end
