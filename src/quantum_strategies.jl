export quantum_strategy

"""
Constructs a strategy matrix given quantum states and measurements. The supported
scenarios include:

[`BlackBox`](@ref) scenarios:

    quantum_strategy(
        Π :: POVM,
        ρ_states :: Vector{<:State};
    ) :: Strategy

For a quantum system the conditional proabilities are constructed as

```math
    P(y|x) = \\text{Tr}[\\Pi_y\\rho_x].
```
"""
function quantum_strategy(
    Π :: POVM,
    ρ_states :: Vector{<:State}
) :: Strategy
    conditionals = measure(Π, ρ_states)
    scenario = BlackBox(length(ρ_states), length(Π))
    Strategy(conditionals, scenario)
end

"""
[`LocalSignaling`](@ref) scenarios:

    quantum_strategy(
        Π :: POVM,
        ρ_states :: Vector{<:State},
        scenario :: LocalSignaling
    ) :: Strategy

For quantum systems the conditional probabilities are construct as

```math
    P(y|x) = \\text{Tr}[\\Pi_y\\rho_x].
```

A `DomainError` is thrown if the provided states and measurements are not compatible
with the specified scenario.
"""
function quantum_strategy(
    Π :: POVM,
    ρ_states :: Vector{<:State},
    scenario :: LocalSignaling
) :: Strategy
    if (size(Π[1]) != (scenario.d, scenario.d)) || (size(ρ_states[1]) != (scenario.d, scenario.d))
        throw(DomainError((Π, ρ_states), "POVM or States are not dimension `d=$(scenario.d)`."))
    end
    conditionals = measure(Π, ρ_states)
    Strategy(conditionals, scenario)
end

"""
[`BipartiteNonSignaling`](@ref) scenarios:

    quantum_strategy(
            ρ_AB :: State,
            Π_Ax :: Vector{<:POVM},
            Π_By :: Vector{<:POVM},
            scenario :: BipartiteNonSignaling;
            atol::Float64 = 1e-7
    )

Constructs a strategy matrix in the generalized representation for the quantum
system with conditional probabilities.

```math
    P(ab|xy) = \\text{Tr}[(\\Pi_a^x\\otimes\\Pi_b^y)\\rho_{AB}]
```

A `DomainError` is thrown if
* The length of each POVM does not match the scenarios number of outputs
* The number of each party's POVMS doesn't match the the number of inputs.
"""
function quantum_strategy(
        ρ_AB :: State,
        Π_Ax :: Vector{<:POVM},
        Π_By :: Vector{<:POVM},
        scenario :: BipartiteNonSignaling;
        atol::Float64 = 1e-7
) :: Strategy
    if scenario.X != length(Π_Ax)
        throw(DomainError(Π_Ax, "`length(Π_Ax) == $(scenario.X)` is not satisfied."))
    elseif scenario.Y != length(Π_By)
        throw(DomainError(Π_By, "`length(Π_By) == $(scenario.Y)` is not satisfied."))
    elseif all(Π -> length(Π) != scenario.A , Π_Ax)
        throw(DomainError(Π_Ax, "Each POVM in Π_Ax must have $(scenario.A) elements"))
    elseif all(Π -> length(Π) != scenario.B, Π_By)
        throw(DomainError(Π_By, "Each POVM in Π_By must have $(scenario.B) elements"))
    end

    strat_mat = zeros(Float64, (scenario.A,scenario.B,scenario.X,scenario.Y))
    for a in 1:scenario.A, b in 1:scenario.B, x in 1:scenario.X, y in 1:scenario.Y
        Π_AB = kron(Π_Ax[x][a].M, Π_By[y][b].M)
        strat_mat[b,a,y,x] = tr(Π_AB * ρ_AB)
    end

    Strategy(reshape(strat_mat, strategy_dims(scenario)), scenario, atol=atol)
end
