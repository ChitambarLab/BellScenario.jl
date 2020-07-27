"""

"""
module Symmetry

using ..BellScenario: PrepareAndMeasure

using ..QBase: QMath

abstract type Group end

function groups(PM :: PrepareAndMeasure)
    input_relabels = QMath.permutation_maps(PM.X)
    output_relabels = QMath.permutation_maps(PM.B)

    Dict(
        "input" => input_relabels,
        "output" => output_relabels
    )
end

end