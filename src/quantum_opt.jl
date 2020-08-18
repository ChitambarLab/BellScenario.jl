using Convex, SCS

export optimize_measurement

"""
    optimize_measurement(
        game :: BellGame,
        ρ_states :: Vector{<:States.AbstractDensityMatrix},
        PM :: PrepareAndMeasure
    )

Perform semi-definite programming to find the optimal quantum measurement which
maximizes the score for the specified set of quantum  states `ρ_states` and `BellGame`.
"""
function optimize_measurement(
    game::BellGame,
    ρ_states::Vector{<:States.AbstractDensityMatrix},
    PM::PrepareAndMeasure,
) :: Dict
    if PM.X != length(ρ_states)
        throw(DomainError(PM, "expexted length of ρ_states is $(PM.d), but got $(length(ρ_states)) instead"))
    end

    if size(ρ_states[1]) != (PM.d,PM.d)
        throw(DomainError(ρ_states, "dimension of ρ_states is not $(PM.d)"))
    end

    norm_game_vector = convert(Vector{Int64}, game)
    norm_bound = norm_game_vector[end]
    norm_game = reshape(norm_game_vector[1:(end-1)], (PM.B-1, PM.X))



    # add povm variables and constraints
    Π_vars = map(i -> HermitianSemidefinite(PM.d), 1:PM.B)
    constraints = (sum(map(Π_b -> real(Π_b), Π_vars)) == Matrix{Float64}(I, PM.d, PM.d))
    constraints += (sum(map(Π_b -> imag(Π_b), Π_vars)) == zeros(Float64, PM.d, PM.d))

    # sum up the state contibutions for each row
    H_b = map(row_id -> sum(norm_game[row_id,:] .* ρ_states), 1:PM.B-1)

    # add the objective
    objective = maximize(real(tr(sum(Π_vars[1:end-1] .* H_b))), constraints)

    # optimize model
    solve!(objective, SCS.Optimizer(verbose=0))

    # parse/return results
    score = objective.optval
    violation = score - norm_bound
    Π_opt = map(Π_b -> Π_b.value, Π_vars)

    Dict(
        "violation" => violation,
        "povm" => Π_opt,
        "game" => game,
        "scenario" => PM,
        "states" => ρ_states
    )
end