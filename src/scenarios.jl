export Scenario, BlackBox, Bipartite, LocalSignaling

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
    )

A black-box signaling scenario involving a transmitting and receiving device.
The transmitter takes `X` inputs and the receiver has `B` outputs.
The transmitter signals to the receiver using `d`-dits of classical or quantum
communication.

![Local Signaling Scenario](../assets/scenario_images/local_signaling_scenario.png)

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

A non-signaling Bell scenario using two devices.

![Bipartite Non-Signaling Scenario](../assets/scenario_images/bipartite_non_signaling_scenario.png)

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
    Bipartite(
        A :: Tuple{Int64, Int64},
        B :: Tuple{Int64, Int64};
        dits :: Int = 1,
        bidirectional :: Bool = false
    ) <: Scenario

A black-box scenario with two devices and possible communication between the
devices. The keyword parameter `dits` describes the number of dits used for
communication and `bidirectional` describes whether communication is simultaneous
between each party.
"""
struct Bipartite <: Scenario
    A :: BlackBox
    B :: BlackBox
    dits :: Int
    bidirectional :: Bool
    Bipartite(
        A::Tuple{Int,Int}, B::Tuple{Int,Int};
        dits::Int=1, bidirectional::Bool=false
    ) = (dits >= 1) ? new(BlackBox(A...), BlackBox(B...), dits, bidirectional) : throw(
            DomainError(dits, "communication `dits` must be ≥ 1.")
        )
end
