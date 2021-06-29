export Scenario, BlackBox, BipartiteSignaling, LocalSignaling

export BipartiteNonSignaling

"""
A `Scenario` is an abstract type parent to all Bell scenarios.
Each child of this abstract type describes a distinct black-box system configuration.
"""
abstract type Scenario end

"""
    BlackBox(num_out :: Int64, num_in :: Int64) <: Scenario

A Bell scenario consisting of a single black-box device with `num_in` inputs and
`num_out` outputs.

![Black-Box Device](../assets/scenario_images/blackbox_device.png)

The black-box device computes the output ``y`` from the input ``x`` by performing
a stochastic map ``S(y|x)``.

### Errors
A `DomainError` is thrown if parameters `num_out` or `num_in` is less than 1.
"""
struct BlackBox <: Scenario
    num_out :: Int64
    num_in :: Int64
    BlackBox(
        num_out::Int64,
        num_in::Int64
    ) = ((num_out >= 1) && (num_in >= 1)) ? new(num_out, num_in) : throw(
        DomainError((num_out, num_in), "inputs must be ≥ 1")
    )
end

"""
    LocalSignaling(
        X :: Int64,
        Y :: Int64,
        d :: Int64,
    ) <: Scenario

A bipartite signaling scenario where information is passed from a transmitter
black-box to a receiver black-box using no more than `d` dits of communication.

![Local Signaling Scenario](../assets/scenario_images/local_signaling_scenario.png)

The transmitter device has `X` inputs and the receiver device has `Y` outputs and
shared randomness is held between the two devices.
When quantum communication is used instead of classical communication no Bell violations
occur.

### Errors
A `DomainError` is thrown if `X`, `Y`, or `d` is less than 1.
"""
struct LocalSignaling <: Scenario
    X :: Int64
    Y :: Int64
    d :: Int64
    LocalSignaling(X::Int64, Y::Int64, d::Int64) = (X ≥ 1 && Y ≥ 1 && d ≥ 1) ? new(X,Y,d) : throw(
        DomainError((X, Y, d), "parameters must be ≥ 1")
    )
end

"""
    BipartiteNonSignaling(
        A :: Int64,
        B :: Int64,
        X :: Int64,
        Y :: Int64
    ) <: Scenario

A bipartite non-signaling scenario where each device receives an input and produces
an output.
Let Alice be the device with `A` outputs and `X` inputs while Bob is the device
with `B` outputs and `Y` inputs.

![Bipartite Non-Signaling Scenario](../assets/scenario_images/bipartite_non_signaling_scenario.png)

Shared randomness is held between Alice and Bob.
When Alice and Bob share quantum entanglement, Bell violations are known to occur.

### Errors
A `DomainError` is thrown if `A`, `B`, `X`, or `Y` is less than 1.
"""
struct BipartiteNonSignaling <: Scenario
    A :: Int64
    B :: Int64
    X :: Int64
    Y :: Int64
    BipartiteNonSignaling(
        A::Int64, B::Int64, X::Int64, Y::Int64
    ) = (
            (A >= 1) && (B >= 1) && (X >= 1) && (Y >= 1)
        ) ? new(A, B, X, Y) : throw(
            DomainError((A, B, X, Y), "inputs must be ≥ 1")
        )
end

"""
    BipartiteSignaling(
        A :: Tuple{Int64, Int64},
        B :: Tuple{Int64, Int64};
        dits :: Int64 = 1,
        bidirectional :: Bool = false
    ) <: Scenario

A bipartite signaling scenario where each device can send a message to the other.

![Bipartite Signaling Scenario](../assets/scenario_images/bipartite_signaling_scenario.png)
"""
struct BipartiteSignaling <: Scenario
    A :: BlackBox
    B :: BlackBox
    dits :: Int64
    bidirectional :: Bool
    BipartiteSignaling(
        A::Tuple{Int64,Int64}, B::Tuple{Int64,Int64};
        dits::Int64=1, bidirectional::Bool=false
    ) = (dits >= 1) ? new(BlackBox(A...), BlackBox(B...), dits, bidirectional) : throw(
            DomainError(dits, "communication `dits` must be ≥ 1.")
        )
end
