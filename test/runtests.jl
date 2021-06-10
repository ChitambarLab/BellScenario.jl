using Test, SafeTestsets

const is_WINDOWS = Sys.iswindows()

println("importing BellScenario.jl")
@time using BellScenario

function _test_runner()
    @time @testset "BellScenario.jl" begin
        @testset "unit tests:" begin
            println("running unit tests.")
            for test in readdir("./test/unit/")
                # run only julia files in test directory
                if occursin(r"^.*\.jl$", test)
                    println("./unit/$test")
                    @time include("./unit/$test")
                end
            end
        end
        @time @testset "integration tests:" begin
            println("running integration tests.")

            println("./test/integration/file_io.jl")
            @time @safetestset "./test/integration/file_io.jl" begin include("integration/file_io.jl") end

            if !is_WINDOWS
                println("./test/integration/LocalPolytope/adjacency_decomposition.jl")
                @time @safetestset "./test/integration/LocalPolytope/adjacency_decomposition.jl" begin include("integration/LocalPolytope/adjacency_decomposition.jl") end
            end

            # for test in readdir("./test/integration/")
            #     # run only julia files in test directory
            #     if occursin(r"^.*\.jl$", test)
            #         println("./integration/$test")
            #         @time include("./integration/$test")
            #     end
            # end
        end
        @testset "regression tests:" begin
            println("running regression tests.")
            for test in readdir("./test/regression/")
                # run only julia files in test directory
                if occursin(r"^.*\.jl$", test)
                    println("./regression/$test")
                    @time include("./regression/$test")
                end
            end
        end
    end
end

# Pkg.test("BellScenario") runs from ./test directory. Development tests from root.
dir = pwd()
if occursin(r".*test$", dir)
    cd(_test_runner, "../")
elseif occursin(r".*BellScenario", dir)
    _test_runner()
else
    error("runtests.jl is currently running from the $(pwd()) directory with contents $(readdir()). runtests.jl must be run from the ./BellScenario.jl or ./BellScenario.jl/test directories.")
end
