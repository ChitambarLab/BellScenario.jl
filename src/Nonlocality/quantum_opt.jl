export optimize_measurement

@doc raw"""
[`LocalSignaling`](@ref) scenario:

    optimize_measurement(
        scenario :: LocalSignaling,
        game :: BellGame,
        ρ_states :: Vector{<:States.AbstractDensityMatrix}
    )

Finds the measurement that optimizes the score of the [`BellGame`](@ref) against
the set of quantum states `ρ_states`.
The optimization is performed with the following semi-definite program:

```math
\begin{aligned}
&\max_{\{\Pi_y\}} \sum_{x,y} G_{x,y} \text{Tr}[\Pi_y \rho_x] \\
&s.t. \quad \sum_y \Pi_y = \mathbb{I}, \quad \Pi_y \geq 0
\end{aligned}
```
"""
function optimize_measurement(
    scenario::LocalSignaling,
    game::BellGame,
    ρ_states::Vector{<:States.AbstractDensityMatrix},
) :: Dict
    if scenario.X != length(ρ_states)
        throw(DomainError(scenario, "expected length of `ρ_states` is $(scenario.X)), but got $(length(ρ_states)) instead"))
    end

    if size(ρ_states[1]) != (scenario.d,scenario.d)
        throw(DomainError(ρ_states, "dimension of `ρ_states` is not $(scenario.d)"))
    end

    norm_game_vector = convert(Vector{Int64}, game)
    norm_bound = norm_game_vector[end]
    norm_game = reshape(norm_game_vector[1:(end-1)], (scenario.Y-1, scenario.X))

    # add povm variables and constraints
    Π_vars = map(i -> HermitianSemidefinite(scenario.d), 1:scenario.Y)
    constraints = (sum(map(Π_y -> real(Π_y), Π_vars)) == Matrix{Float64}(I, scenario.d, scenario.d))
    constraints += (sum(map(Π_y -> imag(Π_y), Π_vars)) == zeros(Float64, scenario.d, scenario.d))

    # sum up the state contibutions for each row
    H_y = map(row_id -> sum(norm_game[row_id,:] .* ρ_states), 1:scenario.Y-1)

    # add the objective
    objective = maximize(real(tr(sum(Π_vars[1:end-1] .* H_y))), constraints)

    # optimize model
    solve!(objective, SCS.Optimizer(verbose=0))

    # parse/return results
    score = objective.optval
    violation = score - norm_bound
    Π_opt = _opt_vars_to_povm(map(Π_y -> Π_y.value, Π_vars))

    Dict(
        "violation" => violation,
        "povm" => Π_opt,
        "game" => game,
        "scenario" => scenario,
        "states" => ρ_states
    )
end

@doc raw"""
[`BipartiteNonSignaling`](@ref) scenario:

    optimize_measurement(
        game :: BellGame,
        scenario :: BipartiteNonSignaling,
        ρ_AB :: States.AbstractDensityMatrix;
        A_POVMs :: Vector{<:Observables.AbstractPOVM},
    ) :: Dict

Find Bob's measurement which optimizes a `BellGame`'s score for the shared quantum
state `ρ_AB` and POVM measurement applied by Alice.
The following semi-definite program optimizes the Bob's POVM:

```math
\begin{aligned}
    &\max_{\{\Pi_b^y\}} \sum_{a,b,x,y} G_{a,b,x,y}\text{Tr}[(\Pi_a^x \otimes \Pi_b^y)\rho_AB] \\
    &s.t. \quad \sum_b \Pi_b^y = \mathbb{I},\quad \Pi_b^y \geq 0 \quad \forall\; y
\end{aligned}
```

    optimize_measurement(
        game :: BellGame,
        scenario :: BipartiteNonSignaling,
        ρ_AB :: States.AbstractDensityMatrix;
        B_POVMs :: Vector{<:Observables.AbstractPOVM},
    ) :: Dict

Find Alice's measurement which optimizes a `BellGame`'s score for the shared quantum
state `ρ_{AB}` and POVM measurement applied by Bob.
The following semi-definite program optimizes the Alice's POVM:

```math
\begin{aligned}
&\max_{\{\Pi_a^x\}} \sum_{a,b,x,y} G_{a,b,x,y}\text{Tr}[(\Pi_a^x \otimes \Pi_b^y)\rho_{AB}] \\
&s.t. \quad \sum_a \Pi_a^x = \mathbb{I},\quad \Pi_a^x \geq 0 \quad \forall \;x
\end{aligned}
```
"""
function optimize_measurement(
    scenario::BipartiteNonSignaling,
    game::BellGame,
    ρ_AB :: States.AbstractDensityMatrix;
    A_POVMs=Vector{Observables.POVM}(undef,0) :: Vector{<:Observables.AbstractPOVM},
    B_POVMs=Vector{Observables.POVM}(undef,0) :: Vector{<:Observables.AbstractPOVM}
) :: Dict
    if length(A_POVMs) > 0
        return _optimize_measurement_B(scenario, game, ρ_AB, A_POVMs)
    elseif length(B_POVMs) > 0
        return _optimize_measurement_A(scenario, game, ρ_AB, B_POVMs)
    else
        throw(DomainError((A_POVMs,B_POVMs), "either `A_POVMs` or `B_POVMs` must be specified."))
    end
end

# """
# Helper function for optimizing Bob's POVM
# """
function _optimize_measurement_B(
    scenario::BipartiteNonSignaling,
    game::BellGame,
    ρ_AB :: States.AbstractDensityMatrix,
    A_POVMs :: Vector{<:Observables.AbstractPOVM}
) :: Dict
    if !(scenario.X == length(A_POVMs))
        throw( DomainError(
            (scenario.X, " != length(A_POVMs)"),
            "A distinct POVM is not specified for input `X` of `BipartiteNonSignaling` scenario"
        ))
    end

    if !(all(i -> scenario.A == length(A_POVMs[i]), 1:scenario.X ))
        throw( DomainError(
            A_POVMs, "Each POVM must have $(scenario.A) outputs."
        ))
    end

    ρ_dim = size(ρ_AB,1)
    Π_A_dim = size(A_POVMs[1][1],1)

    # The dimension of Bob's outputs
    Π_B_dim = (ρ_dim % Π_A_dim == 0) ? Int64(ρ_dim / Π_A_dim) : throw(DomainError(
        ρ_dim, "`size(ρ_Ab)` must equal `size(Π_A)*size(ΠB)`."
    ))

    # Bob's POVM variables to optimize
    B_POVMs = map(y -> map(b -> HermitianSemidefinite(Π_B_dim), 1:scenario.B), 1:scenario.Y)

    # Objective function
    problem = maximize(real(
        sum(x -> sum(y -> sum(a -> sum(b ->
            game[(a-1)*scenario.B+b,(x-1)*scenario.Y+y]*tr(kron(A_POVMs[x][a],B_POVMs[y][b])*ρ_AB),
        1:scenario.B), 1:scenario.A), 1:scenario.Y), 1:scenario.X)
    ))

    # Completeness constraints for Bob's POVM
    for y in 1:scenario.Y
        problem.constraints += sum(real.(B_POVMs[y])) == Matrix{Float64}(I, Π_B_dim, Π_B_dim)
        problem.constraints += sum(imag.(B_POVMs[y])) == zeros(Float64, Π_B_dim, Π_B_dim)
    end

    # optimize model
    solve!(problem, SCS.Optimizer(verbose=0))

    # parse/return results
    score = problem.optval
    violation = score - game.β

    Π_B_opt = map(y -> _opt_vars_to_povm(map(Π_b -> Π_b.value, B_POVMs[y])), 1:scenario.Y)

    Dict(
        "violation" => violation,
        "score" => score,
        "game" => game,
        "scenario" => scenario,
        "state" => ρ_AB,
        "A_POVMs" => A_POVMs,
        "B_POVMs" => Π_B_opt
    )
end

# """
# Helper function for optimizing Alice's POVM
# """
function _optimize_measurement_A(
    scenario::BipartiteNonSignaling,
    game::BellGame,
    ρ_AB :: States.AbstractDensityMatrix,
    B_POVMs :: Vector{<:Observables.AbstractPOVM}
) :: Dict
    if !(scenario.Y == length(B_POVMs))
        throw( DomainError(
            (scenario.Y, " != length(B_POVMs)"),
            "A distinct POVM is not specified for input `Y` of `BipartiteNonSignaling` scenario"
        ))
    end

    if !(all(i -> scenario.B == length(B_POVMs[i]), 1:scenario.Y ))
        throw( DomainError(
            B_POVMs, "Each POVM must have $(scenario.B) outputs."
        ))
    end

    ρ_dim = size(ρ_AB,1)
    Π_B_dim = size(B_POVMs[1][1],1)

    # dimension of Alice's POVM
    Π_A_dim = (ρ_dim % Π_B_dim == 0) ? Int64(ρ_dim / Π_B_dim) : throw(DomainError(
        ρ_dim, "`size(ρ_Ab)` must equal `size(Π_A)*size(ΠB)`."
    ))

    # Alice's POVM variables to optimize
    A_POVMs = map(x -> map(a -> HermitianSemidefinite(Π_A_dim), 1:scenario.A), 1:scenario.X)

    # objective function
    problem = maximize(real(
        sum(x -> sum(y -> sum(a -> sum(b ->
            game[(a-1)*scenario.B+b,(x-1)*scenario.Y+y]*tr(kron(A_POVMs[x][a],B_POVMs[y][b])*ρ_AB),
        1:scenario.B), 1:scenario.A), 1:scenario.Y), 1:scenario.X)
    ))

    # adding completeness constraints
    for x in 1:scenario.X
        problem.constraints += sum(real.(A_POVMs[x])) == Matrix{Float64}(I, Π_A_dim, Π_A_dim)
        problem.constraints += sum(imag.(A_POVMs[x])) == zeros(Float64, Π_A_dim, Π_A_dim)
    end

    # optimize model
    solve!(problem, SCS.Optimizer(verbose=0))

    # parse/return results
    score = problem.optval
    violation = score - game.β

    Π_A_opt = map( x -> _opt_vars_to_povm(map(Π_a -> Π_a.value, A_POVMs[x])), 1:scenario.X)

    Dict(
        "violation" => violation,
        "score" => score,
        "game" => game,
        "scenario" => scenario,
        "state" => ρ_AB,
        "A_POVMs" => Π_A_opt,
        "B_POVMs" => B_POVMs
    )
end

# """
# Helper function for converting POVM optimization variables to a POVM
# """
function _opt_vars_to_povm(povm::Vector{Matrix{Complex{Float64}}}) :: Observables.POVM
    Π = Observables.is_povm(povm) ? Observables.POVM(povm) : throw(DomainError(povm, "not a povm"))

    Π
end
