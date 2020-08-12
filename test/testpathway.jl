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
            @test index(genes, [2, 3]) == [3, 2]
            @test_throws KeyError index(genes, 4)
            @test_throws MethodError index(genes, '3')
            @test_throws KeyError index(genes, [2, 4])
            @test_throws MethodError index(genes, ['2', '3'])
        end
    end

    @testset "Proteins" begin
        ps = [1, 3, 2, 4]
        nin = 1
        nout = 1
        @testset "Constructors" begin
            prtns = Proteins{Int}(ps, Dict(1 => 1, 3 => 2, 2 => 3, 4 => 4), nin, nout)
            prtns_srt = Proteins{Int}(sort(ps), Dict(1 => 1, 2 => 2, 3 => 3, 4 => 4), nin, nout)
            @test Proteins(ps, nin, nout) == prtns
            @test Proteins(ps[begin:begin+nin-1], ps[end-nout+1:end], ps[begin+nin:end-nout]) == prtns
            @test Proteins(Set(ps[begin:end-nout]), Set(ps[begin+nin:end])) == prtns_srt
            @test_throws DomainError Proteins([1, 2, 2, 4], nin, nout)
            @test_throws AssertionError Proteins(ps, 2, 3)
        end
        @testset "Array interface" begin
            prtns = Proteins(ps, nin, nout)
            @test size(prtns) == (4,)
            @test prtns[2] == 3
            @test_throws ErrorException prtns[2] = 5
            @test_throws ErrorException similar(prtns)
            @test length(prtns) == 4
            @test 3 in prtns
            @test 5 ∉ prtns
            @test [p for p in prtns] == [1, 3, 2, 4]
            @test collect(prtns) == [1, 3, 2, 4]
            @test prtns[begin] == 1
            @test prtns[2:end] == [3, 2, 4]
            @test_throws BoundsError prtns[-1]
            @test_throws ErrorException copy(prtns)
        end
        @testset "Other methods" begin
            prtns = Proteins(ps, nin, nout)
            @test prtns === prtns
            @test prtns !== Proteins(ps, nin, nout)
            @test prtns != ps
            @test prtns != "array"
            @test index(prtns, 3) == 2
            @test index(prtns, [2, 3]) == [3, 2]
            @test_throws KeyError index(prtns, 5)
            @test_throws MethodError index(prtns, '3')
            @test_throws KeyError index(prtns, [3, 5])
            @test_throws MethodError index(prtns, ['2', '3'])
            @test input(prtns) == [1,]
            @test output(prtns) == [4,]
            @test activators(prtns) == [1, 3, 2]
            @test products(prtns) == [3, 2, 4]
        end
    end

    @testset "BinaryPhenotype" begin
        prtns = Proteins([1, 3, 2, 4], 1, 1)
        st = BitVector([true, false, true, false])
        pht = BinaryPhenotype(prtns, st)
        @test proteins(pht) === prtns
        @test state(pht, 1) == true
        @test state(pht, 3) == false
        @test_throws KeyError state(pht, 5)
        @test_throws MethodError state(pht, '1')
        @test state(pht, [1, 2, 3]) == BitVector([true, true, false])
        @test_throws KeyError state(pht, [1, 2, 5])
        @test_throws MethodError state(pht, ['1', '2', '3'])
    end

end