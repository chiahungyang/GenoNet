using Test
using GenoNet.PathwayFramework

@testset "PathwayFramework Tests" begin

    @testset "Genes" begin
        gs = [1, 3, 2]
        @testset "Constructors" begin
            @test Genes(gs) == Genes{Int}(gs, Dict(1 => 1, 3 => 2, 2 => 3))
            @test Genes(Set(gs)) == Genes{Int}(sort(gs), Dict(1 => 1, 2 => 2, 3 => 3))
            @test_throws AssertionError Genes([1, 2, 2])
        end
        @testset "Array interface" begin
            gns = Genes(gs)
            @test size(gns) == (3,)
            @test gns[2] == 3
            @test_throws ErrorException gns[2] = 4
            @test_throws ErrorException similar(gns)
            @test length(gns) == 3
            @test 3 in gns
            @test 4 ∉ gns
            @test [g for g in gns] == [1, 3, 2]
            @test collect(gns) == [1, 3, 2]
            @test gns[begin] == 1
            @test gns[2:end] == [3, 2]
            @test_throws BoundsError gns[-1]
            @test_throws ErrorException copy(gns)
        end
        @testset "Other methods" begin
            gns = Genes(gs)
            @test gns === gns
            @test gns !== Genes(gs)
            @test gns != gs
            @test gns != "array"
            @test index(gns, 3) == 2
            @test index(gns, [2, 3]) == [3, 2]
            @test_throws KeyError index(gns, 4)
            @test_throws MethodError index(gns, '3')
            @test_throws KeyError index(gns, [2, 4])
            @test_throws MethodError index(gns, ['2', '3'])
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
            @test_throws AssertionError Proteins([1, 2, 2, 4], nin, nout)
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
        @test isa(pht, AbstractPhenotype)
        @test proteins(pht) === prtns
        @test state(pht, 1) == true
        @test state(pht, 3) == false
        @test_throws KeyError state(pht, 5)
        @test_throws MethodError state(pht, '1')
        @test state(pht, [1, 2, 3]) == BitVector([true, true, false])
        @test_throws KeyError state(pht, [1, 2, 5])
        @test_throws MethodError state(pht, ['1', '2', '3'])
    end

    @testset "BinaryEnv" begin
        prtns = Proteins([1, 3, 2, 4], 1, 1)
        env = BinaryEnv(prtns, [1,], Int[], [4,])
        @testset "Constructors" begin
            @test env.ps === prtns && env.stml == [1,] && env.essnt == Int[] && env.ftl == [4,]
            @test isa(env, AbstractEnv)
            @test_throws AssertionError BinaryEnv(prtns, [1,], Int[], [5,])
            @test_throws MethodError BinaryEnv(prtns, ['1',], Char[], ['4',])
            @test_throws AssertionError BinaryEnv(prtns, [1,], [4,], [4,])
            @test_throws AssertionError BinaryEnv(prtns, [4,], Int[], [4,])
            @test_throws AssertionError BinaryEnv(prtns, [1,], Int[], [1,])
        end
        @testset "Other methods" begin
            @test proteins(env) === prtns
            @test stimulated(env) == [1,]
            @test essential(env) == Int[]
            @test fatal(env) == [4,]
        end
    end

    @testset "DyadicGenotype" begin
        gns = Genes([1, 3, 2])
        prtns = Proteins([1, 3, 2, 4], 1, 1)
        als = Dict(1 => (1 => 3), 2 => (2 => 3), 3 => (3 => 4))
        gt = DyadicGenotype(gns, prtns, als)
        @testset "Constructors" begin
            @test gt.gs === gns && gt.ps === prtns && gt.abstr == [1 => 2, 2 => 4, 3 => 2]
            @test isa(gt, AbstractGenotype)
            @test_throws AssertionError DyadicGenotype(gns, prtns, Dict(1 => (1 => 3), 2 => (2 => 3), 4 => (3 => 4)))
            @test_throws AssertionError DyadicGenotype(gns, prtns, Dict(1 => (1 => 3), 2 => (2 => 3), 3 => (3 => 4), 4 => (3 => 4)))
            @test_throws AssertionError DyadicGenotype(gns, prtns, Dict(1 => (1 => 3), 2 => (2 => 3)))
            @test_throws AssertionError DyadicGenotype(gns, prtns, Dict(1 => (1 => 3), 2 => (2 => 3), 3 => (3 => 5)))
            @test_throws AssertionError DyadicGenotype(gns, prtns, Dict(1 => (1 => 3), 2 => (2 => 3), 3 => (5 => 4)))
            @test_throws AssertionError DyadicGenotype(gns, prtns, Dict(1 => (1 => 3), 2 => (2 => 3), 3 => (3 => 1)))
            @test_throws AssertionError DyadicGenotype(gns, prtns, Dict(1 => (1 => 3), 2 => (4 => 3), 3 => (3 => 4)))
        end
        @testset "Other methods" begin
            @test genes(gt) === gns
            @test proteins(gt) === prtns
            @test allele(gt, 3) == (3 => 4)
            @test_throws KeyError allele(gt, 4)
            @test_throws MethodError allele(gt, '3')
            @test allele(gt, [1, 2, 3]) == [1 => 3, 2 => 3, 3 => 4]
            @test_throws KeyError allele(gt, [1, 2, 4])
            @test_throws MethodError allele(gt, ['1', '2', '3'])
            @test iscompatible(gt, DyadicGenotype(gns, prtns, Dict(1 => (1 => 2), 2 => (3 => 2), 3 => (2 => 4))))
            @test iscompatible(gt, DyadicGenotype(Genes([1, 3, 2]), prtns, als))
            @test iscompatible(gt, DyadicGenotype(gns, Proteins([1, 3, 2, 4], 1, 1), als))
            @test ! iscompatible(gt, DyadicGenotype(Genes([1, 2]), prtns, Dict(1 => (1 => 3), 2 => (2 => 3))))
            @test ! iscompatible(gt, DyadicGenotype(gns, Proteins([1, 2, 3, 4], 1, 1), als))
            env = BinaryEnv(prtns, [1,], Int[], [4,])
            wrongenv = BinaryEnv(Proteins([1, 2, 3, 4], 1, 1), [1,], Int[], [4,])
            @test iscompatible(gt, env)
            @test iscompatible(gt, BinaryEnv(Proteins([1, 3, 2, 4], 1, 1), [1,], Int[], [4,]))
            @test ! iscompatible(gt, wrongenv)
            pht = phenotype(gt, env)
            @test isa(pht, BinaryPhenotype)
            @test pht.ps === prtns && pht.st == BitVector([true, true, false, true])
            @test_throws AssertionError phenotype(gt, wrongenv)
            @test_throws MethodError phenotype(gt, BinaryEnv(Proteins(['1', '4'], 1, 1), ['1',], Char[], ['4',]))
        end
    end

end