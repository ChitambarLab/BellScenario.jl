module Behavior

# Module Description:
#   A behavior is a set of conditional probabilities which describe the relative
#   outcome frequencies for each of unique experimental bell test setting. This
#   module provides utilities for manipulating behaviors for bipartite bell scenarios.
#
#   A behavior can be represented in different subspaces depending on system
#   constraints.

using LinearAlgebra
using QBase: QMath

# """
# Description:
#     Determines the index at which information about the specified outcome is
#     stored.
#
# Inputs:
#     outcome: Tuple, (a,b,x,y) values
#     (α/β)_expt: Tuple, experiment parameters
#     rep: String, behavior representation (generalized,normalized,fixed-direction,no-signaling)
#
# Output:
#     index: Integer, id that indexes into the specified behavior vector
# """
function index(outcome, α_expt, β_expt, rep="generalized")
    (a,b,x,y) = outcome
    (x_in, a_out) = α_expt
    (y_in, b_out) = β_expt

    index = 1
    if rep == "generalized"
        index += (a-1)*b_out*x_in*y_in + (b-1)*x_in*y_in + (x-1)*y_in + y
    elseif rep == "normalized"
        if a == a_out & b == b_out
            throw(ArgumentError("outcome coordinate 'a' or 'b' is not in normalized representation"))
        else
            index += (a-1)*(b_out)*x_in*y_in + (b-1)*x_in*y_in + (x-1)*y_in + y
        end
    elseif rep == "no-signaling"
        if a == a_out | b == b_out
            throw(ArgumentError("outcome coordinate 'a' or 'b' is not in no-signaling representation"))
        elseif (b == 0 & y == 0) & (in(a,1:(a_out-1)) & in(x,1:(x_in)))
            index += (a-1)*x_in + x
        elseif (a == 0 & x == 0) & (in(b,1:(b_out-1)) & in(y,1:(y_in)))
            index += x_in*(a_out - 1) + (b-1)*y_in + y
        elseif (in(a,1:(a_out-1)) & in(x,1:(x_in))) & (in(b,1:(b_out-1)) & in(y,1:(y_in)))
            index += x_in*(a_out - 1) + y_in*(b_out - 1) + (a-1)*(b_out-1)*x_in*y_in + (b-1)*x_in*y_in + (x-1)*y_in + y
        else
            throw(ArgumentError("outcome coordinate ($a,$b,$x,$y) is not valid in no-signaling representation"))
        end
    elseif rep == "fixed-direction"
        if b == b_out
            throw(ArgumentError("outcome coordinate is not in fixed-direction representation"))
        elseif (b == 0 & y == 0) & (in(a,1:(a_out-1)) & in(x,1:(x_in)))
            index += (a-1)*x_in + x
        elseif (in(a,1:a_out) & in(x,1:(x_in))) & (in(b,1:(b_out-1)) & in(y,1:(y_in)))
            index += x_in*(a_out - 1) + (a-1)*(b_out-1)*x_in*y_in + (b-1)*x_in*y_in + (x-1)*y_in + y
        else
            throw(ArgumentError("outcome ($a,$b,$x,$y) is not valid"))
        end
    end

    index
end

# """
# dimension(α_expt,β_expt,rep)
#     Returns the dimension of a behavior vector in the specified representation
#     ("generalized", "normalized", "fixed-direction", "no-signaling").
#
# Inputs:
#     (α/β)_expt: Tuple, bell-test parameters
#     rep: String, behavior representation
#
# Outputs:
#     dim: Integer, behavior dimension
# """
function dimension(α_expt,β_expt,rep)
    (x_num, a_num) = α_expt
    (y_num, b_num) = β_expt
    (num_in, num_out) = α_expt .* β_expt

    dim = 1
    if rep == "generalized"
        dim += num_in * num_out
    elseif rep == "normalized"
        dim += num_in * (num_out - 1)
    elseif rep == "fixed-direction"
        dim += x_num * (a_num - 1) + num_in * a_num * (b_num - 1)
    elseif rep == "no-signaling"
        dim += x_num * (a_num - 1) + y_num * (b_num - 1) + num_in * (a_num - 1) * (b_num - 1)
    end

    dim
end

# """
# Description:
#     Applies tests checking the validity of a behavior in the specified representation.
#
# Input:
#     behavior: Vector
#     (α/β)_expt: experiment parameters
#     rep: String, behavior representation
#
# Output:
#     is_valid: Boolean
# """
function is_valid( behavior, α_expt, β_expt, rep="generalized")
    (x_in, a_out) = α_expt
    (y_in, b_out) = β_expt
    (num_inputs, num_outputs) = α_expt .* β_expt

    is_valid = false

    false_count = 0

    for i in 1:size(behavior)[1]
        if !((0 <= behavior[i]) & (behavior[i] <= 1))
            false_count += 1
        end
    end
    if rep == "generalized"
        # verify normalization
        for x in 1:x_in
            for y in 1:y_in
                norm = 0
                for a in 1:a_out
                    for b in 1:b_out
                        id = index((a,b,x,y),α_expt,β_expt)
                        norm += behavior[id]
                    end
                end

                if !(norm ≈ 1)
                    false_count += 1
                end
            end
        end
    elseif rep == "normalized"
        max_length = Behavior.dimension(α_expt,β_expt,rep) - 1

        valid_checks = map(
            (i) -> begin
                indices = 1 .+ collect(i:num_inputs:max_length)
                norm_sum = sum(behavior[indices])
                (0 <= norm_sum <= 1)
            end,
            1:num_inputs
        )
        if !(all(valid_checks))
            false_count += 1
        end
    elseif rep == "no-signaling"
        #verify normalizability
        for x in 1:x_in
            a_norm = 0
            for a in 1:(a_out-1)
                id = index((a,0,x,0),α_expt,β_expt,"no-signaling")
                a_norm += behavior[id]
            end
            if (1 < a_norm) & !(1 ≈ a_norm)
                false_count += 1
            end
        end

        for y in 1:y_in
            b_norm = 0
            for b in 1:(b_out-1)
                id = index((0,b,0,y),α_expt,β_expt,"no-signaling")
                b_norm += behavior[id]
            end
            if (1 < b_norm) & !(1 ≈ b_norm)
                false_count += 1
            end
        end

        # verify joint construction
        for x in 1:x_in
            for y in 1:y_in
                for a in 1:(a_out-1)
                    for b in 1:(b_out-1)
                        id = index((a,b,x,y),α_expt,β_expt,rep)
                        α_id = index((a,0,x,0),α_expt,β_expt,rep)
                        β_id = index((0,b,0,y),α_expt,β_expt,rep)

                        if !(
                            round(behavior[α_id]*behavior[β_id],sigdigits=3) == round(behavior[id],sigdigits = 3)
                        )
                            false_count += 1
                        end
                    end
                end
            end
        end
    elseif rep == "fixed-direction"
        #verify normalizability
        a_comps = zeros(x_in)
        for x in 1:x_in
            a_norm = 0
            for a in 1:(a_out-1)
                id = index((a,0,x,0),α_expt,β_expt,"fixed-direction")
                a_norm += behavior[id]

            end
            a_comps[x] = 1-a_norm

            if (1 < a_norm) & !(1 ≈ a_norm)
                false_count += 1
            end
        end

        if !all(map((a_comp) -> a_comp >= 0 ,a_comps))
            false_count += 1
        end

        # storing b marginals "separated" from a
        b_margs = fill(-1.0, (x_in,y_in,(b_out-1)))
        for x in 1:x_in
            for y in 1:y_in

                ab_norm = 0
                for a in 1:a_out
                    ax_id = 0
                    ax_val = 0

                    # only test the b_marg if we aren't going to divide by zero
                    test_b_marg = false
                    if a < a_out
                        ax_id = index((a,0,x,0),α_expt, β_expt, "fixed-direction")
                        ax_val = behavior[ax_id]
                        if behavior[ax_id] != 0
                            test_b_marg = true
                        end
                    else
                        if a_comps[x] != 0
                            test_b_marg = true
                        end
                    end

                    for b in 1:(b_out-1)
                        # verifying normalization of joint outputs
                        id = index((a,b,x,y), α_expt, β_expt, "fixed-direction")
                        ab_norm += behavior[id]

                        # verifying that the b_marginal matches the rest
                        if (a < a_out) & test_b_marg
                            if b_margs[x,y,b] == -1
                                b_margs[x,y,b] = (behavior[id]/ax_val)
                            else
                                if !(behavior[id]/ax_val ≈ b_margs[x,y,b])
                                    false_count += 1
                                end
                            end
                        elseif (a == a_out) & test_b_marg
                            if b_margs[x,y,b] == -1
                                b_margs[x,y,b] = behavior[id]/a_comps[x]
                            else
                                if !(behavior[id]/a_comps[x] ≈ b_margs[x,y,b])
                                    false_count += 1
                                end
                            end
                        end

                    end
                end

                if (1 < ab_norm) & !(1 ≈ ab_norm)
                    false_count += 1
                end
            end
        end


    end

    if false_count == 0
        is_valid = true
    end

    is_valid
end

# """
# gen_to_norm_proj(α_expt, β_expt)
#     Constructs the projection matrix from generalized to normalized
#     representation.
#
# Inputs:
#     (α/β)_expt: Tuple, bell-test parameters
#
# Outputs:
#     proj: Matrix that projects gneralized coordinate behaviors to normalized
#        behavior vectors.
# """
function gen_to_norm_proj(α_expt, β_expt)
    (x_num, a_num) = α_expt
    (y_num, b_num) = β_expt
    (num_in, num_out) = α_expt .* β_expt

    gen_dim = num_in * num_out

    proj = zeros(Int64, (1 + num_in*(num_out - 1), gen_dim + 1))

    proj[1,1] = 1

    # fill in 1:1 mappings
    for x in 1:x_num
        for y in 1:y_num
            for a in 1:a_num
                for b in 1:b_num
                    if (a < a_num) | (b < b_num)
                        row = index((a,b,x,y),α_expt,β_expt,"normalized")
                        col = index((a,b,x,y),α_expt,β_expt)
                        proj[row,col] = 1
                    end
                end
            end
        end
    end

    proj
end

# """
# gen_to_fd_proj(α_expt, β_expt)
#     Constructs the projection matrix from generalized to fixed-direction
#     representation.
#
# Inputs:
#     (α/β)_expt: Tuple, bell-test parameters
#
# Outputs:
#     proj: Matrix that projects gneralized coordinate behaviors to fixed-direction
#        behavior vectors.
# """
function gen_to_fd_proj(α_expt, β_expt)
    (x_num, a_num) = α_expt
    (y_num, b_num) = β_expt
    (num_in, num_out) = α_expt .* β_expt

    gen_dim = num_in * num_out
    α_dim = x_num*(a_num - 1)
    β_dim = y_num*(b_num - 1)

    proj = zeros(Int64, (1 + α_dim + x_num*a_num*β_dim, gen_dim + 1))

    proj[1,1] = 1

    # filling in p(a|x) using no-signaling constraint
    for x in 1:x_num
        for a in 1:(a_num - 1)

            row = index((a,0,x,0),α_expt,β_expt,"fixed-direction")

            for b in 1:b_num
                col = index((a,b,x,1),α_expt,β_expt)
                proj[row, col] = 1
            end
        end
    end

    # fill in 1:1 mappings
    for x in 1:x_num
        for y in 1:y_num
            for a in 1:a_num
                for b in 1:(b_num - 1)
                    row = index((a,b,x,y),α_expt,β_expt,"fixed-direction")
                    col = index((a,b,x,y),α_expt,β_expt)
                    proj[row,col] = 1
                end
            end
        end
    end

    proj
end

# """
# gen_to_ns_proj:
#     Returns a projection matrix which maps generalized behaviors to no-signaling
#     behaviors.
#
# Inputs
#     (α/β)_expt: experiment parameters
#
# Output:
#     proj: Matrix which projects generalized to no-signaling representations
# """
function gen_to_ns_proj(α_expt, β_expt)
    (α_num_inputs, α_num_outputs) = α_expt
    (β_num_inputs, β_num_outputs) = β_expt
    (num_inputs, num_outputs) = α_expt .* β_expt

    gen_dim = num_inputs * num_outputs
    α_dim = α_num_inputs*(α_num_outputs - 1)
    β_dim = β_num_inputs*(β_num_outputs - 1)
    joint_dim = α_dim*β_dim

    proj = zeros(Int64, (1 + α_dim + β_dim + joint_dim, 1 + gen_dim))

    proj[1,1] = 1

    # fill in p(a|x)
    for x in 1:α_num_inputs
        for a in 1:(α_num_outputs - 1)

            row = index((a,0,x,0),α_expt,β_expt,"no-signaling")

            for b in 1:β_num_outputs
                col = index((a,b,x,1),α_expt,β_expt)
                proj[row, col] = 1
            end
        end
    end

    # fill in p(b|y)
    for y in 1:β_num_inputs
        for b in 1:(β_num_outputs - 1)

            row = index((0,b,0,y),α_expt,β_expt,"no-signaling")

            for a in 1:(α_num_outputs)
                col = index((a,b,1,y),α_expt,β_expt)
                proj[row, col] = 1
            end
        end
    end

    # fill in 1:1 mappings
    for x in 1:α_num_inputs
        for y in 1:β_num_inputs
            for a in 1:(α_num_outputs - 1)
                for b in 1:(β_num_outputs - 1)
                    row = index((a,b,x,y),α_expt,β_expt,"no-signaling")
                    col = index((a,b,x,y),α_expt,β_expt)
                    proj[row,col] = 1
                end
            end
        end
    end

    proj
end

function norm_to_gen_proj(α_expt, β_expt)
    (x_num, a_num) = α_expt
    (y_num, b_num) = β_expt
    (num_in, num_out) = α_expt .* β_expt

    gen_dim = num_in * num_out

    proj = zeros(Int64, (gen_dim + 1, 1 + num_in*(num_out - 1)))
    proj[1,1] = 1

    for x in 1:x_num
        for y in 1:y_num
            for a in 1:a_num
                for b in 1:b_num
                    # 1:1 mappings
                    if (a < a_num) | (b < b_num)
                        row = index((a,b,x,y),α_expt,β_expt)
                        col = index((a,b,x,y),α_expt,β_expt,"normalized")
                        proj[row,col] = 1
                    else
                        row = index((a_num,b_num,x,y),α_expt,β_expt)
                        proj[row,1] = 1
                        for α in 1:a_num
                            for β in 1:b_num
                                if (α < a_num) | (β < b_num)
                                    col = index((α,β,x,y),α_expt,β_expt,"normalized")
                                    proj[row,col] = -1
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    proj
end

# """
# fd_to_gen_proj(α_expt,β_expt)
#     Constructs the map from fixed-direction to generalized coordinates.
#
# Inputs:
#     (α/β)_expt: Tuple, bell-test parameters
#
# Outputs:
#     proj: Matrix, transforms fixed-directino into generalized coordinates.
# """
function fd_to_gen_proj(α_expt, β_expt)
    (x_num, a_num) = α_expt
    (y_num, b_num) = β_expt
    (num_in, num_out) = α_expt .* β_expt

    gen_dim = num_in * num_out
    α_dim = x_num*(a_num - 1)
    β_dim = y_num*(b_num - 1)

    proj = zeros(Int64, (gen_dim + 1,1 + α_dim + x_num*a_num*β_dim))
    proj[1,1] = 1

    for x in 1:x_num
        for y in 1:y_num

            # 1:1 maps
            for a in 1:a_num
                for b in 1:(b_num - 1)
                    row = index((a,b,x,y),α_expt,β_expt)
                    col = index((a,b,x,y),α_expt,β_expt,"fixed-direction")
                    proj[row,col] = 1
                end
            end

            # p(ak|xy)
            for a in 1:(a_num -1)

                row = index((a,b_num,x,y),α_expt,β_expt)

                α_col = index((a,0,x,0),α_expt,β_expt,"fixed-direction")
                proj[row,α_col] = 1

                for b in 1:(b_num - 1)
                    col = index((a,b,x,y),α_expt,β_expt,"fixed-direction")
                    proj[row, col] = -1
                end
            end

            # populating P(KK|xy)
            row = index((a_num,b_num,x,y),α_expt,β_expt)
            proj[row,1] = 1

            for α in 1:(a_num-1)
                α_col = index((α,0,x,0),α_expt,β_expt,"fixed-direction")
                proj[row,α_col] = -1
            end

            for β in 1:(b_num-1)
                β_col = index((a_num,β,x,y),α_expt,β_expt,"fixed-direction")
                proj[row,β_col] = -1
            end
        end
    end

    proj
end

# """
# ns_to_gen_proj:
#     Returns a projection matrix which maps vectors from the no-signaling representation
#     to the generalized representation.
#
# Inputs
#     (α/β)_expt: experiment parameters
#
# Output:
#     proj: matrix which maps no-signaling to generalized representations
# """
function ns_to_gen_proj(α_expt, β_expt)
    (α_num_inputs, α_num_outputs) = α_expt
    (β_num_inputs, β_num_outputs) = β_expt
    (num_inputs, num_outputs) = α_expt .* β_expt

    gen_dim = num_inputs * num_outputs
    α_dim = α_num_inputs*(α_num_outputs - 1)
    β_dim = β_num_inputs*(β_num_outputs - 1)

    proj = zeros(Int64, (gen_dim + 1,1 + α_dim + β_dim + α_dim*β_dim))
    proj[1,1] = 1

    for x in 1:α_num_inputs
        for y in 1:β_num_inputs
            for a in 1:(α_num_outputs - 1)
                for b in 1:(β_num_outputs - 1)
                    row = index((a,b,x,y),α_expt,β_expt)
                    col = index((a,b,x,y),α_expt,β_expt,"no-signaling")
                    proj[row,col] = 1
                end
            end
        end
    end

    # p(an|xy)
    for x in 1:α_num_inputs
        for y in 1:β_num_inputs
            # p(ak|xy)
            for a in 1:(α_num_outputs -1)

                row = index((a,β_num_outputs,x,y),α_expt,β_expt)

                α_col = index((a,0,x,0),α_expt,β_expt,"no-signaling")
                proj[row,α_col] = 1

                for b in 1:(β_num_outputs - 1)
                    col = index((a,b,x,y),α_expt,β_expt,"no-signaling")
                    proj[row, col] = -1
                end
            end

            # p(kb|xy)
            for b in 1:(β_num_outputs -1)

                row = index((α_num_outputs,b,x,y),α_expt,β_expt)

                β_col = index((0,b,0,y),α_expt,β_expt,"no-signaling")
                proj[row,β_col] = 1

                for a in 1:(α_num_outputs - 1)
                    col = index((a,b,x,y),α_expt,β_expt,"no-signaling")
                    proj[row, col] = -1
                end
            end

            # populating P(KK|xy)
            a = α_num_outputs
            b = β_num_outputs

            row = index((α_num_outputs,b,x,y),α_expt,β_expt)
            proj[row,1] = 1

            for α in 1:(α_num_outputs-1)
                α_col = index((α,0,x,0),α_expt,β_expt,"no-signaling")
                proj[row,α_col] = -1
            end

            for β in 1:(β_num_outputs-1)
                β_col = index((0,β,0,y),α_expt,β_expt,"no-signaling")
                proj[row,β_col] = -1
            end

            for α in 1:(α_num_outputs-1)
                for β in 1:(β_num_outputs-1)
                    col = index((α,β,x,y),α_expt,β_expt,"no-signaling")
                    proj[row,col] = 1
                end
            end
        end
    end

    proj
end

# """
# conditionals(α_expt, β_expt, behavior; rep="generalized")
#
#     Transforms a behavior vector into a conditionals matrix.
#
# Inputs:
#     α/β_expt: Tuple, (num_in, num_out), integer values
#     behavior: Col Vector, a bell-test behavior.
#     rep: String, current behavior representation
#
# Output:
#     conditionals: Matrix, the equivalent conditional probabilities in matrix form.
# """
function conditionals(α_expt, β_expt, behavior; rep="generalized")

    expected_dim = Behavior.dimension(α_expt,β_expt,rep)
    gen_dim = Behavior.dimension(α_expt,β_expt,"generalized")

    (num_in, num_out) = α_expt .* β_expt

    if length(behavior) + 1 == expected_dim
        behavior = reshape(cat([1],behavior,dims=1),(expected_dim,1))
    end

    if !(Behavior.is_valid(behavior, α_expt, β_expt, rep))
        throw(ArgumentError("cannot make conditional probabilities from an invalid behvior"))
    end

    gen_behavior = zeros(gen_dim,1)
    if rep == "generalized"
        gen_behavior = behavior
    elseif rep == "normalized"
        gen_behavior = Behavior.norm_to_gen_proj(α_expt,β_expt)*behavior
    elseif rep == "fixed-direction"
        gen_behavior = Behavior.fd_to_gen_proj(α_expt,β_expt)*behavior
    elseif rep == "no-signaling"
        gen_behavior = Behavior.ns_to_gen_proj(α_expt,β_expt)*behavior
    else
        throw(ArgumentError("invalid behavior representation specified"))
    end

    QMath.Conditionals(Matrix(reshape(gen_behavior[2:end], (num_out,num_in))'))
end

# """
# has_constant(α, β, behavior; rep="normalized")
#
#     Checks if the specified behavior has constant required for performing
#     affine transformations.
#
# Input:
#     α/β: Tuple, specifies inputs and outputs of bipartite bell scenario.
#     behavior: Col vector, may contain the affine constant or not.
#     rep: String, optional, the behavior's representation.
#
# Output:
#     has_constant: Bool, true if the behavior has the correct dimension for
#                   including the constant
# """
function has_constant(α, β, behavior; rep="normalized")
    dim = Behavior.dimension(α, β, rep)
    behavior_dim = length(behavior)

    has_constant = NaN
    if dim == behavior_dim
        has_constant = true
    elseif (dim - 1) == behavior_dim
        has_constant = false
    else
        throw(ArgumentError("invalid behavior dimension"))
    end

    has_constant
end

# """
# add_constants(α, β, behaviors; rep="normalized")
#
#     Appends the constant to each behavior if it does not already exist.
#
# Input:
#     α/β: Tuple, number of inputs and outputs for each party.
#     behaviors: Array, contains behavior column vector
#     rep: String, default "normalized", specifies the representation of the behavior
#
# Output:
#     corrected_behaviors: Array, containsCol Vector, the behavior with the constant
# """
function add_constants(α, β, behaviors; rep="normalized")
    map(
        (behavior) -> begin
            has_constant = Behavior.has_constant(α, β, behavior, rep=rep)

            corrected_behavior = behavior
            if !has_constant
                corrected_behavior = cat([1], behavior, dims=1)
            end

            corrected_behavior
        end,
        behaviors
    )
end

# """
# remove_constants(α, β, behavior; rep="normalized")
#
#     Removes the constant to the behavior if it does not already exist.
#
# Input:
#     α/β: Tuple, number of inputs and outputs for each party.
#     behaviors: Array, contains behavior column vectors.
#     rep: String, default "normalized", specifies the representation of the behavior
#
# Output:
#     corrected_behaviors: Array, contains corrected behavior column vectors.
# """
function remove_constants(α, β, behaviors; rep="normalized")
    map(
        (behavior) -> begin
            has_constant = Behavior.has_constant(α, β, behavior, rep=rep)

            corrected_behavior = behavior
            if has_constant
                corrected_behavior = behavior[2:end]
            end

            corrected_behavior
        end,
        behaviors
    )
end

end
