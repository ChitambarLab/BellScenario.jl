export Scenario, BlackBox, Bipartite, PrepareAndMeasure

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
    PrepareAndMeasure(
        X :: Int,
        B :: Int,
        d :: Int,
    )

A black-box scenario with `X` inputs, `B` outputs, and `d`-dits of communication from the
input device to the output device.
"""
struct PrepareAndMeasure <: Scenario
    X :: Int
    B :: Int
    d :: Int
    PrepareAndMeasure(X::Int, B::Int, d::Int) = (X ≥ 1 && B ≥ 1 && d ≥ 1) ? new(X,B,d) : throw(
        DomainError((X, B, d), "parameters must be ≥ 1")
    )
end
