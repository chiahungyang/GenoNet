using Test
using GenoNet: PathwayFramework

@testset "PathwayFramework Tests" begin

    @testset "Genes" begin
        gs = [1, 3, 2]
        @testset "Constructors" begin
            @test Genes(gs) == Genes{Int}(gs, Dict(1 => 1, 3 => 2, 2 => 3))
            @test Genes(Set(gs)) == Genes{Int}(sort(gs), Dict(1 => 1, 2 => 2, 3 => 3))
            @test_throws DomainError Genes([1, 2, 2])
        end
        @testset "Array interface" begin
            genes = Genes(gs)
            @test size(genes) == (3,)
            @test genes[2] == 3
            @test_throws ErrorException genes[2] = 4
            @test_throws ErrorException similar(genes)
            @test length(genes) == 3
            @test 3 in genes
            @test 4 ∉ genes
            @test [g for g in genes] == [1, 3, 2]
            @test collect(genes) == [1, 3, 2]
            @test genes[begin] == 1
            @test genes[2:end] == [3, 2]
            @test_throws BoundsError genes[-1]
            @test_throws ErrorException copy(genes)
        end
        @testset "Other methods" begin
            genes = Genes(gs)
            @test genes === genes
            @test genes !== Genes(gs)
            @test genes != gs
            @test genes != "array"
            @test index(genes, 3) == 2
            @test_throws KeyError index(genes, 4)
            @test_throws MethodError index(genes, '4')
        end
    end

end