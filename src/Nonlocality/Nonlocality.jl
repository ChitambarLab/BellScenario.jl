"""
*Find the optimal parameters for quantum Bell violations.*

The `Nonlocality.jl` module provides tools to optimize quantum non-locality in
Bell scenarios.

### Module Exports:
* [`optimize_measurement`](@ref): Finds the quantum measurement which violates
    a specified Bell inequality maximally.
"""
module Nonlocality

using ..BellScenario

using QBase, LinearAlgebra

using Convex, SCS

include("./optimize_measurements.jl")

end
