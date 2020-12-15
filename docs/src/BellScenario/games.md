```@meta
CurrentModule = BellScenario
```
# Black-Box Games

```@docs
AbstractGame
Game
BellGame
```

## Conversion Methods

A [`BellGame`](@ref) is a matrix representation of a linear inequality that bounds
a LocalPolytope for some [`Scenario`](@ref). For convenience, conversion methods
are provided to transform BellGames to alternative representations of polytope
facets. There are two supported types which can be converted to/from BellGames:
* Facet - `Vector{Int64}`, a column-major vectorization of the game matrix with  
    the bound placed in the last element of the vector.
* `IEQ` - A facet data structure defined in [XPORTA.jl](https://juliapolyhedra.github.io/XPORTA.jl/dev/).

```@docs
convert(::Type{BellGame},::Vector{Int64},::BlackBox)
convert(::Type{Vector{Int64}}, ::BellGame)
convert(::Type{Vector{BellGame}},::IEQ,::BlackBox)
convert(::Type{IEQ}, bell_games::Vector{BellGame})
```
