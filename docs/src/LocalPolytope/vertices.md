```@meta
CurrentModule = LocalPolytope
```
# Vertices

Vertices are extreme points of the local polytope.
They correspond to deterministic strategies.
The vertex representation of a convex polytope can be computed via the
[Polyhedra.jl](https://github.com/JuliaPolyhedra/Polyhedra.jl) interface.

```@docs
vrep
```

## Vertex Enumeration

The vertices of each local polytope can be enumerated for the specified
Bell [`Scenario`](@ref).

```@docs
vertices
```

## Vertex Counting

Count the number of local polytope vertices for the specified Bell [`Scenario`](@ref).

```@docs
num_vertices
```

## Vertex Dimension

Get the dimension of specified vertex representation.

```@docs
vertex_dims
```
