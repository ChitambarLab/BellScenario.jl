# BellScenario.jl

*Compute Bell inequalities and their quantum violations.*

This software is a general tool for analyzing quantum non-locality in Bell scenarios.

## Featured Modules
  1. [`BellScenario`](@ref): Types and constructors for describing Bell scenarios.
  2. [`LocalPolytope`](@ref): Computes the bounds classical Bell Scenarios.
  3. [`Nonlocality`](@ref): Optimizes quantum states and measurements for non-locality.

# Overview

A Bell scenario describes a system of inter-connected black-box devices.
A black-box is a simple device which takes an input and produces an output, however,
no assumptions are made about the physical process inside the black-box.
Black-boxes and Bell scenarios are completely characterized by their input-output statistics.
The statistics of a Bell scenario depends on the underlying physics used by the black-box devices.

In most familiar cases, a black-box uses *classical* physics.
Some examples of classical black-boxes include cars, slot machines, and software.
Classical physics is deterministic, although, an observer of a classical black-box
may not have access to the parameters required to predict the output.
These parameters are called local hidden variables where local means that events
occurring in one location of space-time do not affect another location.
Essentially, classical black-boxes are deterministic while their local hidden
variables create the perception of probabilistic behavior.
Hence, any probabilistic mixture of deterministic black-box behaviors can, in principle,
be realized for a given Bell scenario.
The bounds of a classical Bell scenario's statistics are regarded as Bell inequalities
and the set of all possible behaviors forms a convex polyhedron known as the *local polytope*.

Bell scenarios are of interest to fundamental physics because not all physical
systems adhere to the constraint of locality.
For example, quantum Bell scenarios are able to realize statistics unattainable
by their corresponding classical systems.
This *non-local* behavior manifests as a strong correlation between distant points
in space-time.
The non-local correlations of quantum systems result in operational advantages in
communications, information security, and data processing applications.

While the success of many future technologies depends on realizing quantum non-locality,
this task is very difficult, even in a theoretical sense.
To witness non-locality, the local bounds must first be derived.
Then, the parameters of the equivalent quantum system must be tuned to violate
the derived Bell inequalities.
Each of these problems are challenging analytically and numerically.
The BellScenario.jl package provides the basic toolkit for finding non-locality
in quantum systems.

# Contents

```@contents
Pages = ["user_guide.md", "BellScenario/overview.md", "LocalPolytope/overview.md", "Nonlocality/overview.md", "development_manual.md",]
Depth = 1
```

# Acknowledgements

Development of BellScenario.jl was made possible by the advisory of Dr. Eric
Chitambar and general support from the Physics Department at the University of
Illinois Urbana-Champaign.
Funding was provided by NSF Award 1914440.

# Citing

To reference this work, see [`CITATION.bib`](https://github.com/ChitambarLab/BellScenario.jl/blob/master/CITATION.bib).

# Licensing

BellScenario.jl is released under the MIT License.
