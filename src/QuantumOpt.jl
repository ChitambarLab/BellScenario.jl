module QuantumOpt

# """
# QuantumOpt:
#
# Module for optimizing quantum measures.
# """

using Convex, SCS, LinearAlgebra

using ..QuantumBehavior

# """
# cvx_prepare_and_measure_max_facet_violation_states(scenario, facet, states):
#
#     Performs semidefinite programming to find the set of quantum measurements
#     which optimally violate the provided facet inequality for the specified
#     quantum states.
#
#   ` This uses cvx optimizer which handles complex-valued hermitian matrices.
#
# Input:
#     num_in/out: Integer, number of inputs/outputs in the scenario
#     facet: Row vector, represents a facet inequality
#     states: Array, contains valid density matrices representing encoded states
#
# Output: Dictionary with keys,
#     violation: Float, the maximal value achieved by quantum strategies. If
#                violation > 0, quantum strategy is nonlocal, violation == 0,
#                strategy lies on facet, violation < 0, strategy lies in polytope
#                volume.
#     povm: array of matrices, the optimal povm.
# """
function cvx_prepare_and_measure_max_facet_violation_states(num_in, num_out, facet, states)
    bound = facet[1]
    dits = size(states[1])[1]

    state_groups = map(i -> sum(facet[(1+1+num_in*(i-1)):(1+num_in+num_in*(i-1))].*states), 1:(num_out - 1))

    # add povm variables and constraints
    povm_variables = map(i -> HermitianSemidefinite(dits), 1:num_out)
    constraints = (sum(map(pv -> real(pv), povm_variables)) == [1 0;0 1])
    constraints += (sum(map(pv -> imag(pv), povm_variables)) == [0 0;0 0])

    # add the objective
    objective = maximize(real(tr(sum(state_groups.*povm_variables[1:(end-1)]))), constraints)

    # optimize model
    solve!(objective, SCS.Optimizer(verbose=0))

    # parse/return results
    violation = bound + objective.optval

    Dict("violation" => violation, "povm" => map(pv -> pv.value, povm_variables))
end

# """
# cvx_prepare_and_measure_max_facet_violation_states2(scenario, facet, states):
#
#     Performs semidefinite programming to find the set of quantum measurements
#     which optimally violate the provided facet inequality for the specified
#     quantum states. This method solves for two povms used.
#
#   ` This uses cvx optimizer which handles complex-valued hermitian matrices.
#
# Input:
#     num_in/out: Integer, number of inputs/outputs in the scenario
#     facet: Row vector, represents a facet inequality
#     states: Array, contains valid density matrices representing encoded states
#
# Output: Dictionary with keys,
#     violation: Float, the maximal value achieved by quantum strategies. If
#                violation > 0, quantum strategy is nonlocal, violation == 0,
#                strategy lies on facet, violation < 0, strategy lies in polytope
#                volume.
#     povms: array of povms, the optimal povm pair.
# """
function cvx_prepare_and_measure_max_facet_violation_states2(α, β, facet, states)
    bound = facet[1]
    dits = size(states[1])[1]

    (α_in, α_out) = α
    (β_in, β_out) = β

    state_groups1 = map(i -> sum(facet[(1+1+α_in*β_in*(i-1)):β_in:(1+α_in*β_in+α_in*β_in*(i-1))].*states), 1:(β_out - 1))
    state_groups2 = map(i -> sum(facet[(1+2+α_in*β_in*(i-1)):β_in:(2+α_in*β_in+α_in*β_in*(i-1))].*states), 1:(β_out - 1))

    povm_variables1 = map(i -> HermitianSemidefinite(dits), 1:β_out)
    povm_variables2 = map(i -> HermitianSemidefinite(dits), 1:β_out)

    constraints = (sum(map(pv -> real(pv), povm_variables1)) == [1 0;0 1])
    constraints += (sum(map(pv -> imag(pv), povm_variables1)) == [0 0;0 0])
    constraints += (sum(map(pv -> real(pv), povm_variables2)) == [1 0;0 1])
    constraints += (sum(map(pv -> imag(pv), povm_variables2)) == [0 0;0 0])

    # add the objective
    objective = maximize(real(tr(sum(state_groups1.*povm_variables1[1:(end-1)] + state_groups2.*povm_variables2[1:(end-1)]))), constraints)

    # optimize model
    solve!(objective, SCS.Optimizer(verbose=0))

    # parse/return results
    violation = bound + objective.optval

    Dict("violation" => violation, "povms" => [map(pv -> pv.value, povm_variables1), map(pv -> pv.value, povm_variables2)])
end

# """
# cvx_prepare_and_measure_max_facet_violation_povm(scenario, facet, povm):
#
#     Performs semidefinite programming to find the set of quantum states
#     which optimally violate the provided facet inequality for the specified
#     quantum measurement.
#
#   ` This uses cvx optimizer which handles complex-valued hermitian matrices.
#
# Input:
#     num_in/out: Integer, number of inputs/outputs in the scenario
#     facet: Row vector, represents a facet inequality
#     povm: Array, contains valid povm to optimize against
#
# Output: Dictionary with keys
#     violation: Float, the maximal value achieved by quantum strategies. If
#                violation > 0, quantum strategy is nonlocal, violation == 0,
#                strategy lies on facet, violation < 0, strategy lies in polytope
#                volume.
#     states: Array of density matrices, the optimal density matrices.
# """
function cvx_prepare_and_measure_max_facet_violation_povm(num_in, num_out, facet, povm)
    bound = facet[1]
    dits = size(povm[1])[1]

    povm_groups = map(i -> sum(facet[(1+i):num_in:end].*povm[1:(end-1)]), 1:num_in)

    # add povm variables and constraints
    ρ_variables = map(i -> HermitianSemidefinite(dits), 1:num_in)
    constraints = map(ρ -> real(tr(ρ)) == 1, ρ_variables)

    # add the objective
    objective = maximize(real(tr(sum(povm_groups.*ρ_variables))), constraints)

    # optimize model
    solve!(objective, SCS.Optimizer(verbose=0))

    # parse/return results
    violation = bound + objective.optval

    Dict("violation" => violation, "states" => map(ρv -> ρv.value, ρ_variables))
end

end
