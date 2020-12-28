using Test

@testset "test/unit/Nonlocality.jl" begin

using BellScenario

@testset "running all tests for module Nonlocality.jl" begin

  println("running Nonlocality.jl unit tests:")
  for test in readdir("./test/unit/Nonlocality/")
      # run only julia files in test directory
      if occursin(r"^.*\.jl$", test)
          println("./unit/Nonlocality/$test")
          @time include("./Nonlocality/$test")
      end
  end

end

end
