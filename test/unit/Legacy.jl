using Test

@testset "test/unit/Legacy.jl" begin

@testset "running all Legacy tests for BellScenario.jl" begin

  println("running Legacy unit tests:")
  for test in readdir("./test/unit/Legacy/")
      # run only julia files in test directory
      if occursin(r"^.*\.jl$", test)
          println("./unit/Legacy/$test")
          @time include("./Legacy/$test")
      end
  end

end

end
