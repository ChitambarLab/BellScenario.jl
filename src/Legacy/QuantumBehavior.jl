module QuantumBehavior

# """
# QuantumBehavior Module
#
# Constructs bell scenario behaviors using the laws of quantum mechanics.
# """

using ..Behavior
using LinearAlgebra

# """
# bipartite_scenario(α_povms, β_povms, ψ)
#
#     Constructs the quantum behavior corresponding to the system described by the
#     provided POVM operators and wavefunction.
#
# Assumption:
#     The number of outputs is uniform across inputs.
#
# Input:
#     α/β_povms: Array, contains the number of POVMs as there are inputs. Each POVM
#                is an array of matrices containing the num-outputs matrices.
#     ψ: Vector, the total wavefunction specified for the bipartite system as a whole.
#
# Output:
#     behavior: Vector, outcome probabilities for all possible events in bell scenario.
#               by default, the no-signaling representation is used.
# """
function bipartite_scenario(α_povms, β_povms, ψ)

    α_expt = (size(α_povms)[1], size(α_povms[1])[1])
    β_expt = (size(β_povms)[1], size(β_povms[1])[1])

    gen_behavior = zeros((Behavior.dimension(α_expt,β_expt,"generalized"),1))
    gen_behavior[1] = 1

    xc = 0
    yc = 0
    ac = 0
    bc = 0
    for x_povm in α_povms
        xc += 1
        for y_povm in β_povms
            yc += 1
            for a in x_povm
                ac += 1
                for b in y_povm
                    bc += 1
                    p_abxy = ψ'*kron(a,b)*ψ

                    index = Behavior.index((ac,bc,xc,yc),α_expt,β_expt,"generalized")

                    gen_behavior[index] = p_abxy
                end
                bc = 0
            end
            ac = 0
        end
        yc = 0
    end
    xc = 0

    behavior = nothing
    if (α_expt[2] < 2) & (β_expt[2] < 2)
        behavior = gen_behavior
    elseif (α_expt[2] < 2) | (β_expt[2] < 2)
        behavior = Behavior.gen_to_norm_proj(α_expt,β_expt)*gen_behavior
    else
        behavior = Behavior.gen_to_ns_proj(α_expt,β_expt)*gen_behavior
    end

    behavior
end

# """
# prepare_and_measure(ρ_set, povms):
#
#     Computes the quantum behavior for a prepare and measure scenario where Alice
#     encodes their input into a quantum state (ρ) and sends it to Bob who then
#     measures the quantum state with the specified POVM.
#
# Input:
#     ρ_set: Array, contains valid density matrices
#     povms: Array, contains a valid povm arrays. One povm for each of Bob's inputs
#
# Output:
#     behavior: Col Vector, a behavior represented in the normalized subspace.
# """
function prepare_and_measure(ρ_set, povms; rep="normalized")

    α_in = length(ρ_set)
    β_in = length(povms)
    num_out = length(povms[1])

    α_expt = (α_in,1)
    β_expt = (β_in,num_out)

    gen_behavior = zeros((Behavior.dimension(α_expt,β_expt,"generalized"),1))
    gen_behavior[1] = 1

    xc = 0
    yc = 0
    bc = 0
    for ρ in ρ_set
        xc += 1
        for povm in povms
            yc += 1
            for b in povm
                bc += 1
                index = Behavior.index((1,bc,xc,yc),α_expt,β_expt,"generalized")
                p_abxy = real(tr(ρ*b))
                gen_behavior[index] = p_abxy
            end
            bc = 0
        end
        yc = 0
    end

    behavior = gen_behavior
    if rep == "normalized"
        norm_proj = Behavior.gen_to_norm_proj(α_expt,β_expt)
        behavior = norm_proj*gen_behavior
    end

    behavior
end

end
