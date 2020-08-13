using Test

@testset "./test/integration/LocalPolytope.jl" begin

  println("running LocalPolytope.jl integration tests:")
  for test in readdir("./test/integration/LocalPolytope/")
      # run only julia files in test directory
      if occursin(r"^.*\.jl$", test)
          println("./integration/LocalPolytope/$test")
          @time include("./LocalPolytope/$test")
      end
  end

end
