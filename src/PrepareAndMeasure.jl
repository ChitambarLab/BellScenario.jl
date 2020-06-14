module PrepareAndMeasure

# """
# PrepareAndMeasure Module
#
#     Generic methods for prepare and measure scenarios
# """

using ..LocalPolytope

# """
# num_zero_rows(strategy):
#
#     Counts the number of empty rows in the strategy matrix.
#
# Input:
#     strategy: Matrix, a valid strategy
#
# Output:
#     num_zero_rows: Int, the number of zero rows.
# """
function num_zero_rows(strategy)
    (rows, cols) = size(strategy)

    length(filter( row -> 0 == sum(strategy[row,:]), 1:rows))
end

# """
# vertices_α_interpolable(strategy1, strategy2):
#
#     Checks if two strategies are linearly interpolable over communication protocols.
#
# Input:
#     strategy1/2: Matrix, valid strategies of the same size.
#
# Output:
#     α_interpolable: Boolean, true if a linear interpolation exists on Alice's side
# """
function vertices_α_interpolable(strategy1, strategy2)
    num_zero1 = num_zero_rows(strategy1)
    num_zero2 = num_zero_rows(strategy2)

    main_strat = (num_zero1 <= num_zero2) ? strategy1 : strategy2
    goal_strat = (num_zero1 <= num_zero2) ? strategy2 : strategy1

    (rows, cols) = size(strategy1)
    strat_perms = LocalPolytope.strategies(cols, rows)

    α_interpolable = 0 < length(filter(perm -> main_strat*perm == goal_strat, strat_perms))

    α_interpolable
end

# """
# vertices_β_interpolable(strategy1, strategy2):
#
#     Checks if two strategies are linearly interpolable over decoding strategiess.
#
# Input:
#     strategy1/2: Matrix, valid strategies of the same size.
#
# Output:
#     β_interpolable: Boolean, true if a linear interpolation exists on Bob's side
# """
function vertices_β_interpolable(strategy1, strategy2)
    num_zero1 = num_zero_rows(strategy1)
    num_zero2 = num_zero_rows(strategy2)

    main_strat = (num_zero1 <= num_zero2) ? strategy1 : strategy2
    goal_strat = (num_zero1 <= num_zero2) ? strategy2 : strategy1

    (rows, cols) = size(strategy1)
    strat_perms = LocalPolytope.strategies(cols, rows)

    β_interpolable = 0 < length(filter(perm -> perm * main_strat == goal_strat, strat_perms))

    β_interpolable
end

# """
# vertices_interpolable(strategy1, strategy2):
#
#     Checks if two strategies are linearly interpolable over encodings OR decodings.
#     If the interpolation is over encodings AND decodings then the interpolation
#     is a non-linear polynomial, the function returns false for this case.
#
# Input:
#     strategy1/2: Matrix, valid strategies of the same size.
#
# Output:
#     interpolable: Boolean, true if a linear interpolation exists for either Alice or Bob.
# """
function vertices_interpolable(strategy1, strategy2)
    (vertices_α_interpolable(strategy1, strategy2) || vertices_β_interpolable(strategy1, strategy2))
end

end # module
