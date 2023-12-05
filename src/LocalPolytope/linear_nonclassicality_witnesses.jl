"""
    linear_nonclassicality_witness(
        vertices :: Vector{Vector{T}} where T <: Real,
        test_behavior :: Vector{<:Real};
        verbose=false :: Bool
    ) :: Vector{Float64}

Obtains a linear inequality ``(\\mathbf{G}^\\star, \\beta^\\star)`` bounding the convex hull of the
set of `vertices` (``\\mathcal{V}``) and is violated by the `test_behavior` ``(\\mathbf{P}\\notin \\text{Conv}(\\mathcal{V}))``.
This task is achieved using the  linear program described by Brunner et al. in Eq. 19 of
[https://arxiv.org/abs/1303.2849](https://arxiv.org/abs/1303.2849). The linear program is


```math
\\begin{align}
    \\min_{(\\mathbf{G}, \\beta)} \\quad & \\langle \\mathbf{G},\\mathbf{P}\\rangle - \\beta & \\notag\\\\
    \\text{s.t.} \\quad  & \\langle \\mathbf{G}, \\mathbf{V} \\rangle - \\beta \\leq 0 & \\quad \\forall \\;\\; \\mathbf{V} \\in \\mathcal{V}  \\notag\\\\
    & \\langle \\mathbf{G}, \\mathbf{P} \\rangle \\leq 1 & \\notag\\\\
\\end{align}
```

The solution to the linear program is a linear nonclassicality witness ``(\\mathbf{G}^\\star, \\beta^\\star)`` becuase the
constraint ``\\lange\\mathbf{G}^\\star, \\mathbf{V} \\rangle - \\beta^\\star \\leq 0`` ensures that no behavior in the polytope
``\\text{Conv}(\\mathcal{V})`` can violate the inequality. Provided that ``\\mathbf{P} \\notin \\text{Conv}(\\mathcal{V})``
technique the output linear inequality witnesses the nonclassicality of the `test_behavior`.

The optimized facet inequality ``(\\mathbf{G}^\\star, \\beta^\\star)`` is returned as a vector ``(G^\\star_{0,0}, \\dots, G^\\star_{Y,X}, -\\beta^\\star)``
where ``G^\\star_{y,x}`` are the elements of ``\\mathbf{G}^\\star``.

!!! note "Supporting Software"
    The linear programming is performed using [HiGHS](https://highs.dev/)
    solver via the [`JuMP`](https://jump.dev/JuMP.jl/stable/)
    interface. Please refer to the source code for more details.

!!! note "Converting Output into Bell Game"
    The linear programming software outputs numerical values that have numerical error. Moreover, the linear inequality is
    scaled such that the classical bound is zero and the `test_behavior` score is one. In order to convert the output
    inequality into a `BellGame`, care must be taken to obtain the correct scaling factor to ensure that elements are integers.

!!! note "Classical Test Behavior"
    If the `test_behavior` ``\\mathbf{P}`` is classical, meaning it satisfies ``\\mathbf{P}\\in\\text{Conv}(\\mathcal{V})``,
    then the optimized linear inequality is the zeros vector.
"""
function linear_nonclassicality_witness(
    vertices :: Vector{Vector{T}} where T <: Real,
    test_behavior :: Vector{<:Real};
    verbose=false :: Bool
) :: Vector{Float64}
    dim_v = length(vertices[1]) + 1

    # initializing model
    model = Model(HiGHS.Optimizer)
    set_attribute(model, MOI.Silent(), !verbose)

    # adding variable for inequality
    @variable(model, s[1:dim_v])

    # adding constraints to model
    for v in vertices
        va = [v..., -1]
        @constraint(model, sum(s.*va) <= 0)
    end

    @constraint(model, c, sum(s.*[test_behavior..., -1]) <= 1)

    # defining the optimization objective
    @objective(model, Max, sum(s.*[test_behavior..., -1]))

    # optimizing
    optimize!(model)

    # return optimized linear inequality in vector form
    return value.(s)
end
