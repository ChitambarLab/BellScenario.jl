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

## Conversion Methods
```@docs
convert(::Type{DeterministicStrategy}, ::Vector{Int64}, ::BlackBox; ::String)
convert(::Type{<:AbstractStrategy}, ::Vector{Float64}, ::BipartiteNonSignaling; ::String)
convert(::Type{Vector{Int64}}, ::DeterministicStrategy; ::String)
```
