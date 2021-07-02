export quantum_strategy

"""
Constructs a strategy matrix given quantum states and measurements. The supported
scenarios include:

[`BlackBox`](@ref) scenarios:

    quantum_strategy(
        Π :: AbstractVector{<:AbstractMatrix},
        ρ_states :: Vector{<:State};
        atol=1e-7::Float64
    ) :: Strategy

For a quantum system the conditional proabilities are constructed as

```math
    P(y|x) = \\text{Tr}[\\Pi_y\\rho_x].
```
"""
function quantum_strategy(
    Π :: AbstractVector{<:AbstractMatrix},
    ρ_states :: Vector{<:AbstractMatrix};
    atol=1e-7::Float64
) :: Strategy
    povm = Π isa POVM ? Π : POVM(Π, atol=atol)
    states = map(ρ -> ρ isa State ? ρ : State(ρ, atol=atol), ρ_states)

    conditionals = measure(povm, states)
    Strategy(conditionals.distribution, atol=atol)
end

"""
[`LocalSignaling`](@ref) scenarios:

    quantum_strategy(
        Π :: AbstractVector{<:AbstractMatrix},
        ρ_states :: Vector{<:AbstractMatrix},
        scenario :: LocalSignaling;
        atol=1e-7::Float64
    ) :: Strategy

For quantum systems the conditional probabilities are construct as

```math
    P(y|x) = \\text{Tr}[\\Pi_y\\rho_x].
```

A `DomainError` is thrown if the provided states and measurements are not compatible
with the specified scenario.
"""
function quantum_strategy(
    Π :: AbstractVector{<:AbstractMatrix},
    ρ_states :: Vector{<:AbstractMatrix},
    scenario :: LocalSignaling;
    atol=1e-7::Float64
) :: Strategy
    if (size(Π[1]) != (scenario.d, scenario.d)) || (size(ρ_states[1]) != (scenario.d, scenario.d))
        throw(DomainError((Π, ρ_states), "POVM or States are not dimension `d=$(scenario.d)`."))
    end

    povm = Π isa POVM ? Π : POVM(Π, atol=atol)
    states = map(ρ -> ρ isa State ? ρ : State(ρ, atol=atol), ρ_states)

    conditionals = measure(povm, states)
    Strategy(conditionals.distribution, scenario, atol=atol)
end

"""
[`BipartiteNonSignaling`](@ref) scenarios:

    quantum_strategy(
            ρ_AB :: AbstractMatrix,
            Π_Ax :: Vector{<:AbstractVector{<:AbstractMatrix}},
            Π_By :: Vector{<:AbstractVector{<:AbstractMatrix}},
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
        ρ_AB :: AbstractMatrix,
        Π_Ax :: Vector{<:AbstractVector{<:AbstractMatrix}},
        Π_By :: Vector{<:AbstractVector{<:AbstractMatrix}},
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

    povm_Ax = map(Π -> Π isa POVM ? Π : POVM(Π, atol=atol), Π_Ax)
    povm_By = map(Π -> Π isa POVM ? Π : POVM(Π, atol=atol), Π_By)
    state = ρ_AB isa State ? ρ_AB : State(ρ_AB, atol=atol)

    strat_mat = zeros(Float64, (scenario.A,scenario.B,scenario.X,scenario.Y))
    for a in 1:scenario.A, b in 1:scenario.B, x in 1:scenario.X, y in 1:scenario.Y
        Π_AB = kron(povm_Ax[x][a].M, povm_By[y][b].M)
        strat_mat[b,a,y,x] = tr(Π_AB * state)
    end

    Strategy(reshape(strat_mat, strategy_dims(scenario)), scenario, atol=atol)
end
