export quantum_strategy

"""
Constructs a strategy matrix given quantum states and measurements. The supported
scenarios include:

`BlackBox` scenarios

    quantum_strategy(
        Π :: Observables.AbstractPOVM,
        ρ_states :: Vector{<:States.AbstractDensityMatrix}
    ) :: Strategy

`PrepareAndMeasure` scenarios

    quantum_strategy(
        Π :: Observables.AbstractPOVM,
        ρ_states :: Vector{<:States.AbstractDensityMatrix},
        PM :: PrepareAndMeasure
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
    PM :: PrepareAndMeasure
) :: Strategy
    dits = PM.d
    if (size(Π[1]) != (dits, dits)) || (size(ρ_states[1]) != (dits, dits))
        throw(DomainError((Π, ρ_states), "POVM or States are not dimension `d=$dits`."))
    end
    conditionals = measurement_probs(Π, ρ_states)
    Strategy(conditionals, PM)
end
