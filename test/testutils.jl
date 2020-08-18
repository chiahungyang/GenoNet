using Test
using GenoNet: Utils

@testset "Utils Tests" begin
    @testset "reachable" begin
        adjlist = [[2, 3], [4,], [4,], Int[]]
        @test Utils.reachable(adjlist, [1,]) == BitVector([1, 1, 1, 1])
        @test Utils.reachable(adjlist, [2,]) == BitVector([0, 1, 0, 1])
        @test Utils.reachable(adjlist, [2, 3]) == BitVector([0, 1, 1, 1])
        @test_throws AssertionError Utils.reachable(adjlist, [5,])
        @test_throws AssertionError Utils.reachable(adjlist, [1, 5])
    end
end
