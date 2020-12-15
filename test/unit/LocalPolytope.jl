using Test, LinearAlgebra

@testset "test/unit/LocalPolytope.jl" begin

using BellScenario

@testset "running all tests for module LocalPolytope.jl" begin

  println("running LocalPolytope.jl unit tests:")
  for test in readdir("./test/unit/LocalPolytope/")
      # run only julia files in test directory
      if occursin(r"^.*\.jl$", test)
          println("./unit/LocalPolytope/$test")
          @time include("./LocalPolytope/$test")
      end
  end

end

end
