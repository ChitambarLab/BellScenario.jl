export Scenario, BlackBox, Bipartite, LocalSignaling

export BipartiteNonSignaling

"""
An abstract type to represent general black-box scenarios.
"""
abstract type Scenario end

"""
    BlackBox(num_in :: Int64, num_out :: Int64) <: Scenario

A black-box scenario considering a single device. A black-box is an uncharacterized
device with a finite number of classical inputs and outputs.

A `DomainError` is throw if parameters `num_in` or `num_out` is less than 1.
"""
struct BlackBox <: Scenario
    num_in :: Int64
    num_out :: Int64
    BlackBox(
        num_in::Int64,
        num_out::Int64
    ) = ((num_in >= 1) && (num_out >= 1)) ? new(num_in, num_out) : throw(
        DomainError((num_in, num_out), "inputs must be ≥ 1")
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

"""
    LocalSignaling(
        X :: Int64,
        B :: Int64,
        d :: Int64,
    )

A black-box signaling scenario involving a transmitting and receiving device.
The transmitter takes `X` inputs and the receiver has `B` outputs.
The transmitter signals to the receiver using `d`-dits of classical or quantum
communication.
"""
struct LocalSignaling <: Scenario
    X :: Int64
    B :: Int64
    d :: Int64
    LocalSignaling(X::Int64, B::Int64, d::Int64) = (X ≥ 1 && B ≥ 1 && d ≥ 1) ? new(X,B,d) : throw(
        DomainError((X, B, d), "parameters must be ≥ 1")
    )
end
