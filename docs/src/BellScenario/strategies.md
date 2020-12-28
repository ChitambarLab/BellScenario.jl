```@meta
CurrentModule = BellScenario
```
# BellScenario.jl - Strategies

```@docs
AbstractStrategy
Strategy
DeterministicStrategy
is_deterministic
strategy_dims
random_strategy
deterministic_strategies
convert(::Type{<:AbstractStrategy}, ::Vector{Float64}, ::BipartiteNonSignaling; ::String)
```

# Quantum Strategies

```@docs
quantum_strategy
```
