"""
The `Nonlocality` submodule provides tools to help find quantum non-locality in
Bell scenarios.

The methods included in this module:
* [`optimize_measurement`](@ref): Finds the quantum measurement which violates
    a specified Bell inequality maximally.
"""
module Nonlocality

using ..BellScenario

using QBase, LinearAlgebra

using Convex, SCS

include("./quantum_opt.jl")

end
