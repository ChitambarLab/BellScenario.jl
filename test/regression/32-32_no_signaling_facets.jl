using Test

@testset "(32)(32) no-signaling regression" begin

using BellScenario: LocalPolytope, ConvexPolytope

# Description:
#   Parses inequalities published at http://www.itp.tuwien.ac.at/~svozil/publ/3-3.html
#   from .txt file. Inequalities are represented in the no-signaling subspace for
#   3 inputs and 2 outputs.
#
# Input: None.
#
# Output: Array, contains 16-dim inequality vectors representing a facet inequality
#
function svozil_pitowsky_inequalities()
    filepath = "./test/regression/files/32-32_svozil_pitowsky_facets.txt"

    (inequalities) = open(filepath) do file
        lines = readlines(file)

        ids = Dict(
            "a1" => 2,
            "a2" => 3,
            "a3" => 4,
            "b1" => 5,
            "b2" => 6,
            "b3" => 7,
            "a1b1" => 8,
            "a1b2" => 9,
            "a1b3" => 10,
            "a2b1" => 11,
            "a2b2" => 12,
            "a2b3" => 13,
            "a3b1" => 14,
            "a3b2" => 15,
            "a3b3" => 16,
        )

        inequalities = []

        for id in 7:2:size(lines)[1]
            line = lines[id]

            lhs_match = match(r"^#\d+:\s(\d)\s\?", line)
            rhs_matches = collect( eachmatch(r"([+-])\s?(\d)?(a\db\d|a\d|b\d|)", line) )

            if size(rhs_matches)[1] > 0
                ineq = zeros(Int64, (1,16))
                ineq[1] = parse(Int64, "-" * lhs_match.captures[1])

                for match in rhs_matches
                    index = ids[match[3]]
                    sign = match[1]
                    multiplier = match.captures[2] === nothing ? "1" : match.captures[2]

                    ineq[index] = parse(Int64, sign * multiplier )
                end

                push!(inequalities, ineq)
            end
        end

        (inequalities)
    end

    inequalities
end

@testset "verifying polytope" begin
    dir = "./test/regression/files/"

    vertices = LocalPolytope.vertices((3,2),(3,2))

    @test size(vertices[1]) == (15,1)

    constraints = ConvexPolytope.facet_constraints(vertices, dir=dir)
    @test constraints[2] == []
    @test constraints[3] == []

    inequalities = constraints[1]
    @test size(inequalities) == (684,)
    @test size(inequalities[1]) == (1,16)
    @test issetequal(inequalities, svozil_pitowsky_inequalities())
end

end
