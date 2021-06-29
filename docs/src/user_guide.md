# User Guide

## Quickstart

```julia
julia> using Pkg; Pkg.add("BellScenario")
```

```@example tutorial
using BellScenario
```

## CHSH Scenario

The CHSH scenario is a [`BipartiteNonSignaling`](@ref) scenario where Alice and Bob
each have a black-box with binary inputs and outputs.

![Classical CHSH Scenario](../assets/scenario_images/classical_chsh_scenario.png)

This scenario is significant because it is the simplest Bell scenario in which
quantum nonlocality can be observed.
We will use BellScenario.jl to compute the CH Bell inequality and optimize quantum
measurements to violate the CH inequality.
First, create a CHSH `Scenario` to specify the black-box arrangement in the figure
above.

```@example tutorial
# (num_out_A, num_out_B, num_in_A, num_in_B)
chsh_scenario = BipartiteNonSignaling(2,2,2,2)
```

Bell inequalities bound the set of local (classical) correlations.
The local set is a convex polytope referred to as the *local polytope* and the
facets of the local polytope are Bell inequalities.
The standard method of computing Bell inequalities is to first compute the local
polytope vertices, then apply a polytope transformation algorithm to compute the
Bell inequalities.

The BellScenario.jl package provides the [`LocalPolytope`](@ref) module to compute
Bell inequalities.
The first step is to enumerate the vertices for the CHSH scenario.

```@example tutorial
chsh_vertices = LocalPolytope.vertices(chsh_scenario)
```

Then, the Bell inequalities can computed using the [`LocalPolytope.facets`](@ref)
function.

```@example tutorial
chsh_facets = LocalPolytope.facets(chsh_vertices)["facets"]
```

We'll take ``15^{th}`` facet as it represents the CH inequality

```math
- P_A(1|2) - P_B(1|1) + P(11|11) - P(11|12) + P(11|21) + P(11|22) \leq 0.
```

In fact, this inequality is equivalent to the more celebrated CHSH inequality.
The difference is that the CH inequality is expressed in terms of probabilities
whereas the CHSH inequality is expressed in terms of bipartite correlators.

```@example tutorial
ch_inequality = chsh_facets[15]
```

Now that we have computed a Bell inequality, we can find a quantum violation using
the [`Nonlocality`](@ref) module.
This will require the use of the [QBase.jl](https://github.com/ChitambarLab/QBase.jl) package.
In this example, we will fix Alice's measurement and the quantum state shared
between Alice and Bob.

```@example tutorial
using QBase

# maximally entangled state
ρ_AB = State([1 0 0 1;0 0 0 0;0 0 0 0;1 0 0 1]/2)

# Alice's measurement
A_POVMs = [
    POVM([ [1 0;0 0], [0 0;0 1] ]),
    POVM([ [1 1;1 1]/2, [1 -1;-1 1]/2 ])
]
```

Then, we convert the `ch_inequality` into a general representation of a [`BellGame`](@ref).

```@example tutorial
ch_game = convert(BellGame, ch_inequality, chsh_scenario)
```

Finally, we optimize Bob's measurement with respect to the fixed state and measurements.

```@example tutorial
opt_dict = Nonlocality.optimize_measurement(
    chsh_scenario, ch_game, ρ_AB, A_POVMs=A_POVMs
)
```

We see that the inequality is violated for the optimized measurement and states.

```@example tutorial
opt_dict["violation"]
```
