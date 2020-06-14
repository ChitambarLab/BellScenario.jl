using Test

@testset "running tests for module BellComm.jl" begin

    println("running BellComm.jl unit tests:")
    for test in readdir("./test/BellComm/unit/")
        # run only julia files in test directory
        if occursin(r"^.*\.jl$", test)
            println("./BellComm/unit/$test")
            include("./unit/$test")
        end
    end

    println("running BellComm.jl integration tests:")
    for test in readdir("./test/BellComm/integration/")
        if occursin(r"^.*\.jl$", test)
            println("./BellComm/integration/$test")
            include("./integration/$test")
        end
    end
end
