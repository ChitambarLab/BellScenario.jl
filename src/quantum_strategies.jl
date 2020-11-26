export quantum_strategy

"""
Constructs a strategy matrix given quantum states and measurements. The supported
scenarios include:

`BlackBox` scenarios

    quantum_strategy(
        Π :: Observables.AbstractPOVM,
        ρ_states :: Vector{<:States.AbstractDensityMatrix}
    ) :: Strategy

`LocalSignaling` scenarios

    quantum_strategy(
        Π :: Observables.AbstractPOVM,
        ρ_states :: Vector{<:States.AbstractDensityMatrix},
        scenario :: LocalSignaling
    ) :: Strategy

A `DomainError` is thrown if the provided states and measurements are not compatible
with the specified scenario.
"""
function quantum_strategy(
    Π :: Observables.AbstractPOVM,
    ρ_states :: Vector{<:States.AbstractDensityMatrix}
) :: Strategy
    conditionals = measurement_probs(Π, ρ_states)
    scenario = BlackBox(length(ρ_states), length(Π))
    Strategy(conditionals, scenario)
end

function quantum_strategy(
    Π :: Observables.AbstractPOVM,
    ρ_states :: Vector{<:States.AbstractDensityMatrix},
    scenario :: LocalSignaling
) :: Strategy
    if (size(Π[1]) != (scenario.d, scenario.d)) || (size(ρ_states[1]) != (scenario.d, scenario.d))
        throw(DomainError((Π, ρ_states), "POVM or States are not dimension `d=$(scenario.d)`."))
    end
    conditionals = measurement_probs(Π, ρ_states)
    Strategy(conditionals, scenario)
end
