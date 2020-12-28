using XPORTA: IEQ

export AbstractGame, Game, BellGame

"""
An `AbstractGame` is the abstract type that is parent to type representing a cost
function for strategies.
A `Game` is a linear inequality dual to `Strategy` matrices, it is represented by
a matrix ``G`` containing scalar coefficients and a bound ``\\beta``.
The bound of the linear inequality is typically a maximum score attainable by Bell
scenario using a certain set of resources.
A `Game` ``G`` is played by a `Strategy` ``S`` to achieve a score computed as,

```math
\\beta \\geq \\langle G, S\\rangle = \\sum_{x,y} G_{y,x}S(y|x).
```

The game is "won" if the strategy scores greater than ``\\beta``, that is, the above
inequality is violated.
Since ``\\beta`` is the maximum score for a set of resources, the game is won only
if the tested strategy used a set of resources of greater operational value than
than considered when computing the bound ``\\beta``.

Conveniently, `Strategy` matrices have normalized columns and non-negative elements.
This means that any game inequality can be converted into a form where game matrix
``G`` has positive elements and ``\\beta`` designates a positive upper bound.
"""
abstract type AbstractGame{T} <: AbstractMatrix{T} end
Base.size(G::AbstractGame) = size(G.game)
Base.getindex(G::AbstractGame, I::Vararg{Int64,2}) = getindex(G.game, I...)
Base.setindex!(G::AbstractGame, v, I::Vararg{Int64,2}) = (G.game[I...] = v)

"""
    Game(game::Matrix{T}, β::Real) <: AbstractGame{T}

A `Game` is represented by a `Matrix` of coefficients and a scalar bound `β`.

Type parameter `T` can be either `Int64` or `Float64`.
"""
struct Game{T} <: AbstractGame{T}
    game :: Matrix{T}
    β :: Real
    scenario :: Scenario
    Game(game::Matrix{T}, β::Real) where T = new{T}(game, β)
end

"""
    BellGame(game::Matrix{Int64}, β::Int64, scenario::Scenario)

A `BellGame` represents a Bell inequality or tight facet of the local polytope.
Since the vertices of the local polytope are deterministic strategies with 0,1 elements,
the linear inequalities describing facets of the local polytope have rational coefficients.
Therefore, if a inequality tightly bounds the local polytope, it can be represented by
a game with integer coefficients. 
"""
struct BellGame <: AbstractGame{Int64}
    game :: Matrix{Int64}
    β :: Int64
    scenario :: Scenario
    BellGame(game::Matrix{Int64}, β::Int64) = new(game, β)
end
