using Test

@testset "./test/integration/Legacy.jl" begin

  println("running Legacy integration tests:")
  for test in readdir("./test/integration/Legacy/")
      # run only julia files in test directory
      if occursin(r"^.*\.jl$", test)
          println("./integration/Legacy/$test")
          @time include("./Legacy/$test")
      end
  end

end
