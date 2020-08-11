using Test

@testset "./src/file_io.jl" begin

using BellScenario

@testset "pretty_print_txt()" begin
    test_dir = "./test/integration/files/"
    filename = test_dir * "pretty_print_txt_test"
    filename_txt = filename * ".txt"

    if isfile(filename_txt)
        rm(filename_txt)
    end

    bell_games = BellGame.([[1 0 0;1 0 0;0 0 0],[1 0 0;0 1 0;0 0 1]], [1,2])
    match_string = "(1)  1 ≥ \n\t[1, 0, 0]\n\t[1, 0, 0]\n\t[0, 0, 0]\n\n(2)  2 ≥ \n\t[1, 0, 0]\n\t[0, 1, 0]\n\t[0, 0, 1]\n\n"

    @testset "no .txt extension in `filename`" begin
        @test !isfile(filename_txt)

        pretty_print_txt(bell_games, filename)

        @test isfile(filename_txt)

        file_string = read(filename_txt, String)

        @test file_string == match_string

        rm(filename_txt)
    end

    @testset ".txt extension in `filename`" begin
        @test !isfile(filename_txt)

        pretty_print_txt(bell_games, filename_txt)

        @test isfile(filename_txt)

        file_string = read(filename_txt, String)

        @test file_string == match_string

        rm(filename_txt)
    end
end

end
