"""

"""
module Symmetry

using ..BellScenario: LocalSignaling

using ..QBase: QMath

abstract type Group end

function groups(scenario :: LocalSignaling)
    input_relabels = QMath.permutation_matrices(scenario.X)
    output_relabels = QMath.permutation_matrices(scenario.B)

    Dict(
        "input" => input_relabels,
        "output" => output_relabels
    )
end

end
