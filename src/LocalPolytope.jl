module LocalPolytope

# """
# Module Description:
#   This module houses functions used to compute local polytopes for bipartite
#   bell scenarios via the method of vertex (deterministic behavior) enumeration.
# """

using QBase: QMath
# using IterTools

# """
# strategies(num_inputs, num_outputs):
#   Constructs the set of unique strategies which maps M unique inputs to K
#   unique outputs.
#
# Input:
#   num_inputs: Integer, the number of distinct inputs, M
#   num_outputs: Integer, the number of distinct outputs, K
#
# Output:
#   strategies: Array of (M x K) matrices which deterministically and uniquely
#               map an input to an output.
# """
function strategies(num_inputs, num_outputs)

    inputs = QMath.id(num_inputs)
    outputs = QMath.id(num_outputs)

    strategies = []
    if num_outputs > 1
        # loop through the number of ways to select 10 elements from the output set
        # Incrementing a base (num_outputs) number with num_inputs digit will achieve
        # all such combinations.
        for i in 0:(num_outputs^num_inputs - 1)
            base_n_array = digits(i, base = num_outputs, pad = num_inputs )

            # use the base (num_outputs) number to construct the strategy matrix
            output_index_array = base_n_array .+ 1

            push!(strategies,outputs[:,output_index_array])
        end
    else
        strat = fill(1,(1,num_inputs))
        push!(strategies, strat)
    end

    strategies
end

# """
# communication_protocols(num_inputs,num_outputs):
#   Partitions the m input set into k non-empty output sets. For m >= k,
#   the number of partitions are exactly described by stirling's number
#   of the second kind. For m <= k, the number of protocols is given by (k choose m).
#
# Inputs:
#   num_inputs: Integer, the number of distinct inputs, m.
#   num_outputs: Integer, the number of distinct outputs, k.
#
# Output:
#   protocols: Array, where elements are deterministic maps from input states
#              to output states.
# """
function communication_protocols(num_inputs,num_outputs)
    outputs = QMath.id(num_outputs)
    protocols = []

    if num_inputs < num_outputs

        num_base = max(num_outputs,2)

        inputs = QMath.id(num_inputs)

        min_n_array = collect(0:(num_inputs-1))
        min_num = QMath.base_n_val(min_n_array, num_base)

        base_n_array = reverse(digits(min_num, base = num_base, pad = num_inputs))
        id_array = base_n_array .+ 1

        push!(protocols, outputs[:,id_array])
    else
        num_base = max(num_outputs,2)

        max_n_array = cat(
            0:(num_outputs-1),
            fill(num_outputs-1,(num_inputs - num_outputs)),
            dims=1
        )
        max_num = QMath.base_n_val(max_n_array, num_base)

        min_n_array = cat(
            fill(0, (num_inputs - num_outputs)),
            0:(num_outputs - 1),
            dims = 1
        )
        min_num = QMath.base_n_val(min_n_array, num_base)

        partitions = []
        for i in min_num:max_num
            base_n_array = digits(i, base = num_base, pad = num_inputs)
            # subsets are required to be non-empty thus every digit should be represented
            if issetequal( unique(base_n_array), 0:(num_outputs - 1))
                partition = fill([],num_outputs)
                for j in 1:size(base_n_array)[1]
                    digit = base_n_array[j]
                    partition[digit + 1] = cat(partition[digit + 1], j, dims = 1)
                end

                # ordering of sets does not matter
                # [[1],[2],[3,4]] is same as [[2],[1],[4,3]]
                dupe_count = 0
                for set in partitions
                    if issetequal(partition,set[2])
                        dupe_count += 1
                    end
                end

                # only consider partition sets that are unique
                if dupe_count == 0
                    push!(partitions, (base_n_array, partition) )
                end
            end
        end

        # construct protocols from partitions
        for partition in partitions
            output_index_array = reverse(partition[1]) .+ 1
            push!(protocols, outputs[:,output_index_array])
        end
    end

    protocols
end

# """
# input_set(num_inputs):
#   Constructs a set of orthonormal vectors which represent distinct inputs to
#   the bell test.
#
# Input:
#   num_inputs: Integer, the number of distinct inputs
#
# Output:
#   input_set: Array, contains input vectors of dimension num_inputs
# """
function input_set(num_inputs)
    input_set = []
    for i in 1:num_inputs
        input = zeros(Int64, (num_inputs,1))
        input[i] = 1
        push!(input_set, input)
    end

    input_set
end

# """
# is_no_signaling_duplicate(strategy, num_inputs, num_rx_inputs):
#   The set of strategies that have an rx_input contains strategies that which
#   do not use the rx_input to map the output. These strategies yield redundant
#   vertices in a quantity equivalent to Stirling's number of the 2nd kind to
#   those generated by the no-signaling constraint. This function checks if a
#   strategy is a duplicate.
#
# Inputs:
#   strategy: Matrix, map from (num_inputs)⊗(num_rx_inputs) -> num_outputs
#   num_inputs: Integer, number of inputs supplied to expt (M)
#   num_rx_inputs: Integer, number of inputs supplied by communication (2^bits)
#
# Outputs:
#   is_duplicate: Boolean, true if strategy is a no-signaling duplicate
# """
function is_no_signaling_duplicate(strategy, num_inputs, num_rx_inputs)
    # count the times that changing the rx bit changes the output
    is_duplicate = true
    for i in 0:(num_inputs - 1)
        output = strategy[:,(i*num_rx_inputs + 1)]
        for j in 2:num_rx_inputs
            if output != strategy[:,(i*num_rx_inputs + j)]
                is_duplicate = false
                break
            end
        end
        if !is_duplicate
            break
        end
    end

    is_duplicate
end

# """
# compute_no_signaling_behavior(alice_strategy, bob_strategy; subspace="no-signaling"):
#   Given deterministic strategies for alice and bob, computes the no-signaling
#   behavior of the joint system.
#
# Inputs:
#   alice_strategy: Tuple, (strategy map, set of input, input dimension, output dimension)
#   bob_strategy: Tuple, (strategy map, set of input, input dimension, output dimension)
#   subspace: String, "no-signaling" is optimal representation, "normalization" is general representation
#
# Output:
#   behavior: Array, behavior vector in the specified subspace
# """
function compute_no_signaling_behavior(alice_strategy, bob_strategy; subspace="no-signaling")

    (α, α_in_set, α_in, α_out) = alice_strategy
    (β, β_in_set, β_in, β_out) = bob_strategy

    if subspace == "no-signaling"
        # Alice is assumed to be indpendent of Bob
        α_behavior = zeros(Int64, α_in*(α_out - 1))
        for x in α_in_set
            a = α[1:(end-1),:]*x
            α_behavior += kron(a,x)
        end

        β_behavior = zeros(Int64, β_in*(β_out - 1))
        for y in β_in_set
            b = β[1:(end-1),:]*y
            β_behavior += kron(b,y)
        end

        joint_behavior = zeros(Int64, α_in*(α_out - 1)*β_in*(β_out - 1))
        for x in α_in_set
            for y in β_in_set
                a = α[1:(end-1),:]*x
                b = β[1:(end-1),:]*y

                joint_behavior += kron(a, kron(b, kron(x, y)))
            end
        end

        behavior = cat(α_behavior, β_behavior, joint_behavior, dims = 1)
    elseif subspace == "fixed-direction"
        α_behavior = zeros(Int64, α_in*(α_out - 1))
        αβ_behavior = zeros(Int64, α_in*α_out*β_in*(β_out - 1))
        α_complement = zeros(Int64, α_in)
        for x in α_in_set

            # Alice output
            a = α*x

            α_behavior += kron(a[1:(end-1),:], x)

            for y in β_in_set

                b = β[1:(end-1),:]*y
                αβ_behavior += kron(a,kron(b,kron(x,y)))
            end
        end

        behavior = cat(
            α_behavior,
            αβ_behavior,
            dims = 1
        )

    elseif subspace == "normalized"
        behavior_dim = (α_in * α_out)*(β_in * β_out) - (α_in * β_in)
        behavior = zeros(Int64, behavior_dim)
        joint_strategy = kron(α, β)

        for x in α_in_set
            for y in β_in_set
                xy = kron(x,y)

                ab = joint_strategy[1:(end-1),:]*xy

                behavior += kron(ab,xy)
            end
        end
    else
        behavior_dim = (α_in * α_out)*(β_in * β_out)
        behavior = zeros(Int64, behavior_dim)
        joint_strategy = kron(α, β)
        for x in α_in_set
            for y in β_in_set
                xy = kron(x,y)

                ab = joint_strategy*xy

                behavior += kron(ab,xy)
            end
        end
    end

    behavior
end

# """
# compute_signaling_behavior(alice_strategy, bob_strategy; direction="fixed"):
#   Computes the behavior associated wiht the provided strategies and communication
#   protocols. Supports fixed-direction (Assume Alice -> Bob) and one-way (Alice
#   -> Bob or Alice <- Bob).
#
# Inputs
#   alice_strategy: Tuple, (strategy map, input set, num inputs, num outputs, communication protocol map)
#   bob_strategy: Tuple, (strategy map, input set, num inputs, num outputs, communication protocol map)
#   direction: String, "fixed", "<-", or "->"
#
# Output:
#   behavior: Vector, behavior represented in either the fixed-direction
#             represention or normalization representation
# """
function compute_signaling_behavior(alice_strategy, bob_strategy; direction="fixed", subspace="fixed-direction")
    (α, α_in_set, α_in, α_out, ρ) = alice_strategy
    (β, β_in_set, β_in, β_out, σ) = bob_strategy

    if (direction == "fixed") & (subspace != "fixed-direction")
        direction = "->"
    end

    behavior_dim = 0
    if subspace == "normalized"
        behavior_dim = (α_in * α_out)*(β_in * β_out) - α_in*β_in
    elseif subspace == "generalized"
        behavior_dim = (α_in*α_out)*(β_in*β_out)
    elseif subspace == "fixed-direction"
        behavior_dim = α_in*(α_out - 1) + α_in*(α_out)*β_in*(β_out - 1)
    end

    if direction == "fixed"
        # assume alice only signals to bob
        behavior_dim = α_in*(α_out - 1) + α_in*(α_out)*β_in*(β_out - 1)

        α_behavior = zeros(Int64, α_in*(α_out - 1))
        αβ_behavior = zeros(Int64,  α_in*(α_out)*β_in*(β_out - 1))

        behavior = zeros(Int64, behavior_dim)
        for x in α_in_set

            # Alice output
            a = α*x

            α_behavior += kron(a[1:(end-1),:], x)
            for y in β_in_set

                α_tx = ρ*x
                β_input = kron(y, α_tx)

                b = β[1:(end-1),:]*β_input
                αβ_behavior += kron(a, kron(b, kron(x,y)))
            end
        end

        behavior = cat(
            α_behavior,
            αβ_behavior,
            dims = 1
        )
    elseif direction == "->"
        # TODO: compute generalized behavior and normalized behavior based off subpsace
        behavior = zeros(Int64, behavior_dim)
        for x in α_in_set
            for y in β_in_set

                rx_input = ρ*x
                β_input = kron(y, rx_input)


                joint_input = kron(x, β_input)
                joint_strategy = kron(α,β)
                if subspace == "normalized"
                    ab = joint_strategy[1:(end-1),:]*joint_input
                    behavior += kron(ab,kron(x,y))
                elseif subspace == "generalized"
                    ab = joint_strategy*joint_input
                    behavior += kron(ab,kron(x,y))
                end
            end
        end
    elseif direction == "<-"

        behavior = zeros(Int64, behavior_dim)
        for x in α_in_set
            for y in β_in_set

                rx_input = σ*y
                α_input = kron(x, rx_input)


                joint_input = kron(α_input, y)
                joint_strategy = kron(α,β)

                if subspace == "normalized"
                    ab = joint_strategy[1:(end-1),:]*joint_input
                    behavior += kron(ab,kron(x,y))
                elseif subspace == "generalized"
                    ab = joint_strategy*joint_input
                    behavior += kron(ab,kron(x,y))
                end
            end
        end
    end

    behavior
end

# """
# vertices(alice_expt, bob_expt; bits=0, fixed_direction=true):
#   Enumerates the vertices of the local polytope corresponding to the supplied
#   bipartite bell experiment. Supports no-signaling (bits = 0), fixed direction
#   signaling, and one-way signaling scenarios.
#
# Input:
#   alice_expt:      Tuple, (num_inputs, num_outputs)
#   bob_expt:        Tuple, (num_inputs, num_outputs)
#   bits:            Integer, the number of bits used in communication (default = 0)
#   dits:            Integer, the number of distinct communication states (default = 0)
#   fixed_direction: Boolean, true if fixed-direction signaling is used (default = true)
#
# Output:
#   vertices: Array containing local polytope vertices represented in the optimal
#             subspace.
# """
function vertices(alice_expt, bob_expt; bits=0, dits=0, fixed_direction=true, rep="")

    (α_in, α_out) = alice_expt
    (β_in, β_out) = bob_expt

    comm_dits = 0
    if bits > 0
        comm_dits = 2^bits
    elseif dits > 0
        comm_dits = dits
    end

    α_in_set = input_set(α_in)
    β_in_set = input_set(β_in)

    α_set = strategies(α_in, α_out)
    β_set = strategies( ((comm_dits == 0) ? 1 : comm_dits)*β_in, β_out)

    α2_set = []
    β2_set = []
    ρ_set = []
    σ_set = []

    if !in(rep,("normalized","generalized","fixed-direction","no-signaling"))
        if comm_dits == 0
            rep = "no-signaling"
        elseif (comm_dits > 0) & fixed_direction
            rep = "fixed-direction"
        elseif (comm_dits > 0) & !fixed_direction
            rep = "normalized"
        else
            rep = "generalized"
        end
    end

    if (α_out < 2) & (β_out < 2)
        rep = "generalized"
    elseif (α_out < 2) | (β_out < 2)
        if rep != "generalized"
            rep = "normalized"
        end
    end

    if comm_dits > 0
        β2_set = strategies( β_in, β_out )
        ρ_set = communication_protocols(α_in, comm_dits)
        if !(fixed_direction)
            α2_set = strategies( (comm_dits)*α_in, α_out )
            σ_set = communication_protocols(β_in, comm_dits)
        end
    end

    vertices = []

    if comm_dits == 0
        for α in α_set
            for β in β_set
                vertex = compute_no_signaling_behavior(
                    (α, α_in_set, α_in, α_out),
                    (β, β_in_set, β_in, β_out),
                    subspace = rep
                )
                push!(vertices,vertex)
            end
        end
    elseif fixed_direction
        for α in α_set
            for β2 in β2_set
                vertex = compute_no_signaling_behavior(
                    (α, α_in_set, α_in, α_out),
                    (β2, β_in_set, β_in, β_out),
                    subspace = rep
                )
                push!(vertices, vertex)
            end

            for β in β_set
                if !( is_no_signaling_duplicate(β, β_in, comm_dits) ) & !(α_in < 2)
                    for ρ in ρ_set
                        vertex = compute_signaling_behavior(
                            (α, α_in_set, α_in, α_out, ρ),
                            (β, β_in_set, β_in, β_out, nothing),
                            direction = "fixed",
                            subspace = rep
                        )
                        push!(vertices, vertex)
                    end
                end
            end
        end
    else
        for α in α_set
            for β2 in β2_set
                no_signaling_vertex = compute_no_signaling_behavior(
                    (α, α_in_set, α_in, α_out),
                    (β2, β_in_set, β_in, β_out),
                    subspace = rep
                )
                push!(vertices, no_signaling_vertex)
            end
        end

        for β in β_set
            if !( is_no_signaling_duplicate(β, β_in, comm_dits) ) & !(α_in < 2)
                for α in α_set
                    for ρ in ρ_set
                        ab_vertex = compute_signaling_behavior(
                            (α, α_in_set, α_in, α_out, ρ),
                            (β, β_in_set, β_in, β_out, nothing),
                            direction = "->",
                            subspace = rep
                        )
                        push!(vertices, ab_vertex)
                    end
                end
            end
        end

        for α in α2_set
            if !( is_no_signaling_duplicate(α, α_in, comm_dits) ) & !(β_in < 2)
                for β in β2_set
                    for σ in σ_set
                        ba_vertex = compute_signaling_behavior(
                            (α, α_in_set, α_in, α_out, nothing),
                            (β, β_in_set, β_in, β_out, σ),
                            direction = "<-",
                            subspace = rep
                        )

                        push!(vertices, ba_vertex)
                    end
                end
            end
        end
    end

    vertices
end

# """
# num_prepare_and_measure_vertices(N, d)
#
#     Counts the polytope vertices for prepare and measure scenarios where the
#     number of inputs matches that of the outputs.
#
# Input:
#     N: Integer, number of inputs/outputs
#     d: Integer, number of communication dits
#
# Output:
#     Integer, number of vertices
# """
function num_prepare_and_measure_vertices(N, d)
    sum(map(i -> QMath.stirling2(N, i)*binomial(N, i)*factorial(i), 1:d))
end

# """
# num_inhomogeneous_prepare_and_measure_vertices(X, B, d)
#
#     Counts the polytope vertices for prepare and measure scenarios where the
#     number of inputs is different than the number of outputs.
#
# Input:
#     X: Integer, number of inputs
#     B: Integer, number of outputs
#     d: Integer, number of communication dits
#
# Output:
#     Integer, number of vertices
# """
function num_inhomogeneous_prepare_and_measure_vertices(X, B, d)
    sum(map(i -> QMath.stirling2(X, i)*binomial(B, i)*factorial(i), 1:d))
end

# """
# strategy_to_behavior(strategy):
#
#     Converts a strategy matrix to its corresponding behavior representation.
#
# Input:
#     strategy: Matrix, dimensions do not matter, it should represent a valid strategy
#
# Output:
#     behavior: Col Vector, corresponding behavior representation for the strategy
# """
function strategy_to_behavior(strategy; rep="generalized", constant=true)

    if rep == "normalized"
        strategy = strategy[1:(end-1),:]
    end

    (num_out, num_in) = size(strategy)

    behavior_dim = num_out*num_in

    behavior = zeros(behavior_dim, 1) + collect(Iterators.flatten(map(i->strategy[i,:], 1:num_out)))

    if constant
        behavior = cat([1],behavior,dims=1)
    end

    behavior
end

# """
# behavior_to_strategy(num_in, num_out, gen_behavior):
#
#     converts a behavior in the generalized representation into it's correpsonding
#     strategy matrix.
#
# Inputs:
#     num_in/out: Integer, bell scenario parameters (dimensions of strategies)
#     gen_behavior: Col Vector, a valid behavior in the generalized representation
#
# Output:
#     strategy: Matrix, the strategy isomorphic to the provided behavior.
# """
function behavior_to_strategy(num_in, num_out, gen_behavior)
    behavior = []
    if length(gen_behavior) == num_out*num_in
        behavior = gen_behavior
    elseif length(gen_behavior) - 1 == num_out*num_in
        behavior = gen_behavior[2:end]
    else
        throw(DomainError(gen_behavior, "./src/BellComm/LocalPolytope.jl: generalized behavior does not have proper dimensions"))
    end

    strategy = reshape(behavior, num_in, num_out)'

    strategy
end

"""
    facet_to_matrix(num_in, num_out, gen_facet)

transforms the vector form of a facet to the minimal matrix with non-negative
elements.
"""
function facet_to_matrix(num_in, num_out, gen_facet)

    bound = -1*gen_facet[1]

    matrix = behavior_to_strategy(num_in, num_out, gen_facet[2:end])

    new_matrix = zeros(num_out,num_in)
    for col_id in 1:num_in
        col = matrix[:,col_id]

        col_min = min(col...)

        new_col = col
        if col_min != 0
            new_col = -1*col_min .+ col

            bound = bound + -1*col_min
        end

        new_matrix[:,col_id] = new_col
    end

    (bound, new_matrix)
end

end
