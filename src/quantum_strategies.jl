export quantum_strategy

"""
Constructs a strategy matrix given quantum states and measurements. The supported
scenarios include:

`BlackBox` scenarios

    quantum_strategy(
        Π :: POVM,
        ρ_states :: Vector{<:State}
    ) :: Strategy

`LocalSignaling` scenarios

    quantum_strategy(
        Π :: POVM,
        ρ_states :: Vector{<:State},
        scenario :: LocalSignaling
    ) :: Strategy

A `DomainError` is thrown if the provided states and measurements are not compatible
with the specified scenario.
"""
function quantum_strategy(
    Π :: POVM,
    ρ_states :: Vector{<:State}
) :: Strategy
    conditionals = measure(Π, ρ_states)
    scenario = BlackBox(length(ρ_states), length(Π))
    Strategy(conditionals, scenario)
end

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
