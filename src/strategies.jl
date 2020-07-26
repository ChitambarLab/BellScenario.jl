export Strategy, strategy_dims

"""
A strategy matrix describing the statistical behavior of a black-box.

    Strategy(conditionals :: Matrix{<:Real}) <: AbstractMatrix{Float64}

By default, the constructor creates a strategy for a 'BlackBox' scenario. However,
a `Scenario` can also be passed to the `Strategy` constructor.

    Strategy(conditionals :: Matrix{<:Real}, scenario :: Scenario)

A `DomainError` is thrown if the provided `Scenario` does not match the dimension
of the `conditionals` matrix.
"""
struct Strategy <: AbstractMatrix{Float64}
    conditionals ::  QMath.Conditionals
    Base.size(S::Strategy) = size(S.conditionals)
    Base.getindex(S::Strategy, I::Vararg{Int,2}) = getindex(S.conditionals, I...)
    Base.setindex!(S::Strategy, v, I::Vararg{Int,2}) = (S.conditionals[I...] = v)

    scenario :: Scenario
    Strategy(
        conditionals :: Matrix{<:Real},
        scenario::Scenario
    ) =  size(conditionals) == strategy_dims(scenario) ? new(QMath.Conditionals(conditionals), scenario) : throw(
        DomainError(conditionals, "conditionals are not the correct dimension for the specified bell scenario.")
    )
    Strategy( conditionals :: Matrix{<:Real} ) = new(
        QMath.Conditionals(conditionals),
        BlackBox(reverse(size(conditionals))...)
    )
end

"""
Two strategies can be multiplied together resulting in a new strategy, e.g. `S1*S2 = S3`.

    *(S1::Strategy, S2::Strategy) :: Strategy

The chained strategies can describe a new black-box scenario, therefore, `scenario`
parameter can also be passed to the multiplication operator.

    *(S1::Strategy, S2::Strategy, scenario::Scenario) :: Strategy
"""
*(S1::Strategy, S2::Strategy) = Strategy(S1.conditionals*S2.conditionals)
*(S1::Strategy, S2::Strategy, scenario::Scenario) = Strategy(S1.conditionals*S2.conditionals, scenario)

"""
    strategy_dims( scenario :: Scenario ) :: Tuple{Int, Int}

Returns the dimensions of the `Conditionals` matrix describing a `Strategy` for
the `Scenario` at hand. Each `Scenario`, can place unique constraints on the matrix
dimensions, therfore, separate methods are called for each concrete `Scenario`.
"""
function strategy_dims(scenario::BlackBox) :: Tuple{Int, Int}
    (scenario.num_out, scenario.num_in)
end

function strategy_dims(scenario::Bipartite) :: Tuple{Int, Int}
    A = scenario.A
    B = scenario.B

    ( A.num_out * B.num_out, A.num_in * B.num_in )
end
