```@meta
CurrentModule = BellScenario
```
# BellScenario.jl - Strategies

```@docs
AbstractStrategy
Strategy
strategy_dims
random_strategy
```

## Deterministic Strategies

```@docs
DeterministicStrategy
is_deterministic
deterministic_strategies
```

## Quantum Strategies

```@docs
quantum_strategy
```

## Behaviors

The conditional probabilities of a Bell scenario can be represented by data structures
other than a `Strategy` matrix.
Most notably, the conditional probabilities can be organized into a vector referred
to as a *behavior* in the literature.
A behavior is a column-major vectorization of a strategy matrix e.g. `S[:]`.
The dimension of the behavior vector (or strategy matrix) can be reduced to an
equivalent subspace by using the normalization and non-signaling constraints to
removed redundant parameters of the conditional probability distribution.

Subspaces used within BellScenario.jl include:
* `"generalized"` - All conditional probabilities ``P(y|x)`` for the Bell scenario are present.
* `"normalized"` - The normalization constraint on each column is used to remove the last row of the strategy matrix.
* `"non-signaling"` - The non-signaling constraints are applied to reduce the `"normalized"` subspace further.

## Conversion Methods

Strategy matrices and behavior vectors are isomorphic representations.
Therefore, a behavior can be converted into a strategy and vice versa.

Currently, behaviors are loosely represented by `Vector{Int64}` and `Vector{Float64}` types
and a `rep` argument is used to specify the subspace.
This may change in future updates of BellScenario.jl.

```@docs
convert(::Type{DeterministicStrategy}, ::Vector{Int64}, ::BlackBox; ::String)
convert(::Type{<:AbstractStrategy}, ::Vector{Float64}, ::BipartiteNonSignaling; ::String)
convert(::Type{Vector{Int64}}, ::DeterministicStrategy; ::String)
```
