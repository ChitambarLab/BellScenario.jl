module DichotomicLocalPolytope

# """
# Module: DichotomicLocalPolytope
#
# Hard-coded computation of local polytope vertices for the CHSH scenario with
# and without auxiliary communication.
# """

using LinearAlgebra
# """
# bipartite_no_signaling_strategies(;subspace="p16"):
#   Returns the set of strategies for dichotomic, bi-partite bell test. The
#   joint strategies for Alice and Bob are represented as the kronecker product
#   of their individual strategies α(x)⊗β(y).
#
# Input: subspace, string, the behavior subspace for which results are returned
#        defaults to "p16", the most general and high dimensional bipartite space.
#
# Output: An array of joint strategies which map input x⊗y → a⊗b
# """
function bipartite_no_signaling_strategies(;subspace="p16")
    # Assumption: Alice and Bob have dichotomic setups
    s = dichotomic_strategies()
    n = size(s)[1]

    no_signaling_strategies = [] # initialize empty array

    # Alice and Bob each choose a strategy from the dichotomic strategies.
    for α in s
        for β in s
            if subspace == "p16"
                push!( no_signaling_strategies, kron(α, β) )
            elseif subspace == "p12"
                push!( no_signaling_strategies, kron(α, β)[1:3,:] )
            end
        end
    end

    no_signaling_strategies
end

# """
# deterministic_no_signaling_behaviors(;subspace="p16"):
#   Computes the deterministic behaviors (vertices) of the local polytope for
#   a bi-partitie, dichotomic setup for which non-signaling is applied.
#
# Input: subspace, string, the behavior subspace for which results are returned
#        defaults to "p16", the most general and high dimensional bipartite space.
#
# Output: Column Vectors, which represent deterministic behaviors in the 16-dim
#   basis described by { P(ab|xy) } ∀ a,b,x,y.
# """
function deterministic_no_signaling_behaviors(;subspace="p16")
    # Assumptions: 2-parties with dichotomic setups
    if subspace == "p16"
        s = bipartite_no_signaling_strategies()
    elseif subspace == "p12"
        s = bipartite_no_signaling_strategies(subspace="p12")
    end

    n = size(s)[1]

    behaviors = []
    # looping through distinct strategies
    for i in 1:n
        # initializing an empty probabilistic behavior
        if subspace == "p16"
            p = zeros(Int64,16)
        elseif subspace == "p12"
            p = zeros(Int64,12)
        end

        # columns of identity matrix represent xy ∈ {00, 01, 10, 11}
        identity = Matrix(1I,(size(s[i])[2], size(s[i])[2]))

        # Strategy columns represent outputs (a⊗b) for given input (x⊗y). We
        # compute the 16-dim behavior by (a⊗b)⊗(x⊗y) → P(ab|xy)
        for j in 1:size(s[i])[2]
            p = p + kron(s[i][:,j],identity[:,j])
        end

        push!(behaviors,p)
    end

    behaviors
end

# """
# Description:
#   Computes the vertices of the local polytope which bounds the system where
#   Alice sends 1-bit of communication to Bob.
#
# Input:
#   subspace: String (optional), the behavior subspace in which we represent the
#             deterministic behaviors. default, 16-dim
#
# Output:
#   an array of deterministic vertices represented in the specified subspace
# """
function deterministic_fixed_direction_signaling_behaviors(;subspace="p16")
    # Assumption: signaling goes from Alice to Bob

    ds = dichotomic_strategies()

    ts = two_input_strategies()

    binary_inputs = [
        [1;0], # zero
        [0;1]  # one
    ]

    deterministic_behaviors = []
    for α in ds
        for β in ts

            if subspace == "p16"
                behavior = zeros(Int64,16)
            elseif subspace == "p12"
                behavior = zeros(Int64,12)
            elseif subspace == "p10"
                behavior = zeros(Int64,10)
            end

            for x in binary_inputs

                for y in binary_inputs

                    xy = kron(x,y)
                    input = kron(xy,x)

                    if subspace == "p16"
                        ab = kron(α,β)*input
                        behavior += kron(ab,xy)
                    elseif subspace == "p12"
                        ab = kron(α,β)[1:3,:]*input
                        behavior += kron(ab,xy)
                    elseif subspace == "p10"
                        ab = kron(α,β)[[1,3],:]*input

                        behavior += cat(zeros(Int64,2),kron(ab,xy),dims=1)
                    end
                end
                if subspace == "p10"
                    a = α[1,:]'*x # α gets cast as a 1-D column vector
                    behavior += cat(kron(a,x),zeros(Int64,8),dims=1)
                end

            end
            push!(deterministic_behaviors, behavior)
        end
    end

    deterministic_behaviors
end

# """
# deterministic_one_way_signaling_behaviors(;subspace="p16"):
#   Computes the deterministic behaviors (vertices) of the local polytope by
#   executing a finite set of deterministic strategies for all inputs. In this
#   computation it is assumed that 1-bit of information is passed either from
#   Alice to Bob or from Bob to Alice. The case is also considered where there
#   is non-signaling between Alice and Bob.
#
# Input: subspace, string, the behavior subspace for which results are returned
#        defaults to "p16", the most general and high dimensional bipartite space.
#
# Output: An array of p12 behaviors
# """
function deterministic_one_way_signaling_behaviors(;subspace="p16")
    ds = dichotomic_strategies()

    # filter out strategies that do not depend on second input
    ts = filter(f -> (
            (f*[1;0;0;0] != f*[0;1;0;0]) | (f*[0;0;1;0] != f*[0;0;0;1])
        ),
        two_input_strategies()
    )

    binary_inputs = [[1;0],[0;1]]
    deterministic_behaviors = deterministic_no_signaling_behaviors(subspace=subspace)
    # case for when Alice receives 1-bit from Bob
    for α in ts
        for β in ds
            if subspace == "p16"
                strategy_behavior = zeros(Int64,16)
            elseif subspace == "p12"
                strategy_behavior = zeros(Int64,12)
            end

            for x in binary_inputs
                for y in binary_inputs
                    xy = kron(x,y)
                    input = kron(xy,y)

                    if subspace == "p16"
                        ab = kron(α,β)*input
                    elseif subspace == "p12"
                        ab = kron(α,β)[1:3,:]*input
                    end

                    strategy_behavior += kron(ab,xy)
                end
            end
            push!(deterministic_behaviors,strategy_behavior)
        end
    end

    # case for when Bob receives 1-bit from Alice
    for α in ds
        for β in ts
            if subspace == "p16"
                strategy_behavior = zeros(Int64,16)
            elseif subspace == "p12"
                strategy_behavior = zeros(Int64,12)
            end

            for x in binary_inputs
                for y in binary_inputs
                    xy = kron(x,y)
                    input = kron(xy,x)

                    if subspace == "p16"
                        ab = kron(α,β)*input
                    elseif subspace == "p12"
                        ab = kron(α,β)[1:3,:]*input
                    end

                    strategy_behavior += kron(ab,xy)
                end
            end
            push!(deterministic_behaviors,strategy_behavior)
        end
    end

    deterministic_behaviors
end

# """
# dichotomic_strategies():
#   Returns the set of deterministic strategies that map a binary
#   input to a binary output.
#
# Input: None.
#
# Output: An array of 2x2 matrices that represent all possible strategies.
# """
function dichotomic_strategies()
    I2 = [1 0; 0 1]
    PC0 = [1 1; 0 0]
    PC1 = [0 0; 1 1]
    P2 = [0 1; 1 0]

    [I2, PC0, PC1, P2]
end

# """
# two_input_strategies():
#   Returns the set of deterministic strategies that map two binary inputs to
#   a single binary output.
#
# Input: None.
#
# Output: An array of 2x4 matrices that represent all possible strategies
# """
function two_input_strategies()
    s = dichotomic_strategies()
    n = size(s)[1]

    two_input_strategies = []
    for i in 1:n
        for j in 1:n
            push!(two_input_strategies, cat(s[i],s[j],dims=2))
        end
    end

     two_input_strategies
end

end
