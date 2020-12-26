# BellScenario.jl

*Compute Bell inequalities and their quantum violations.*

This software is a general tool for analyzing quantum non-locality in Bell scenarios.

## Featured Modules
  1. [`BellScenario`](@ref): Types and constructors for describing Bell scenarios.
  2. [`LocalPolytope`](@ref): Computes the bounds classical Bell Scenarios.
  3. [`Nonlocality`](@ref): Optimizes quantum states and measurements for non-locality.

# Overview

A Bell scenario describes a black-box system characterized by its input-output
statistics. A black-box is a device that takes an input and produces an output, however,
no other assumptions are made about the physical process inside the device.

When classical physics are applied within a black-box system, the resulting statistics
can be achieved using deterministic devices and randomness.
The statistics form a convex polytope which regarded as the *local polytope*.

Non-locality is a characteristic of quantum theory that allows strong correlations
to be made between distant black-boxes without any communication.
Classical devices are unable to replicate the statistic of quantum non-locality
which leads to advantages in communications, data processing, and information security.

# Acknowledgements

Development of BellScenario.jl was made possible by the advisory of Dr. Eric
Chitambar and general support from the Physics Department at the University of
Illinois Urbana-Champaign.
Funding was provided by NSF Award 1914440.

# Citing

To reference this work, see [`CITATION.bib`](https://github.com/ChitambarLab/BellScenario.jl/blob/master/CITATION.bib).

# Licensing

BellScenario.jl is released under the MIT License.

# Contents

```@contents
Pages = ["user_guide.md", "BellScenario/overview.md", "LocalPolytope/overview.md", "Nonlocality/overview.md", "development_manual.md",]
Depth = 1
```
