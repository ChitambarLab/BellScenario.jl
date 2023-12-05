"""
    linear_nonclassicality_witness(
        vertices :: Vector{Vector{T}} where T <: Real,
        test_behavior :: Vector{<:Real};
        verbose=false :: Bool
    ) :: Vector{Float64}

Implements a linear program that obtains a linear inequality that witnesses the nonclassciality of
the provided `test_behavior` (see https://arxiv.org/abs/1303.2849 Eq. 19). A complete set of
enumerated vertices are required as input. The optimized facet inequality is returned in vector
form with the classical bound being appended to the end.
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
