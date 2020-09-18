using Test
using GenoNet: Utils
using SparseArrays

@testset "Utils Tests" begin

    @testset "reachable" begin
        adjlist = [[2, 3], [4,], [4,], Int[]]
        @test Utils.reachable(adjlist, [1,]) == BitVector([1, 1, 1, 1])
        @test Utils.reachable(adjlist, [2,]) == BitVector([0, 1, 0, 1])
        @test Utils.reachable(adjlist, [2, 3]) == BitVector([0, 1, 1, 1])
        @test_throws AssertionError Utils.reachable(adjlist, [5,])
        @test_throws AssertionError Utils.reachable(adjlist, [1, 5])
    end

    @testset "leadingeigenvec" begin
        mat = spzeros(2, 2)
        mat[1, 1] = 10
        mat[2, 2] = 2
        eigvec = Utils.leadingeigenvec(mat)
        @test isapprox(abs.(eigvec), [1, 0], atol=1e-8)
    end

end
