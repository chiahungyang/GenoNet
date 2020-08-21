using Test
using GenoNet.PathwayFramework
using GenoNet.PopGenModel
using Random

@testset "PopGenModel" begin
    gns = Genes([1, 3, 2])
    prtns = Proteins([1, 3, 2, 4], 1, 1)
    als_viable = Dict(1 => (1 => 3), 2 => (2 => 3), 3 => (2 => 4))
    gt_viable = DyadicGenotype(gns, prtns, als_viable)
    als_inviable = Dict(1 => (1 => 3), 2 => (2 => 3), 3 => (3 => 4))
    gt_inviable = DyadicGenotype(gns, prtns, als_inviable)
    env = BinaryEnv(prtns, [1,], Int[], [4,])

    @testset "BinaryViability" begin
        @test isviable(phenotype(gt_viable, env), env, BinaryViability())
        @test isviable(gt_viable, env, BinaryViability())
        @test ! isviable(phenotype(gt_inviable, env), env, BinaryViability())
        @test ! isviable(gt_inviable, env, BinaryViability())
        wrongenv = BinaryEnv(Proteins([1, 2, 3, 4], 1, 1), [1,], Int[], [4,])
        @test_throws AssertionError isviable(gt_viable, wrongenv, BinaryViability())
    end

    @testset "IdenticalReproductivity" begin
        popl = [gt_viable, gt_inviable, gt_inviable]
        offsp = similar(popl)
        @test begin
            Random.seed!(12345)
            reproduce!(offsp, popl, IdenticalReproductivity())
            offsp == [gt_inviable, gt_inviable, gt_viable]
        end
        @test begin
            Random.seed!(12345)
            reproduce(popl, 3, IdenticalReproductivity()) == [gt_inviable, gt_inviable, gt_viable]
        end
    end

    @testset "IndependentMutation" begin
        @test begin
            Random.seed!(12345)
            gt_mut = mutate(gt_inviable, IndependentMutation(0.5))
            gt_mut == DyadicGenotype(gns, prtns, Dict(1 => (1 => 3), 2 => (1 => 2), 3 => (3 => 4)))
        end
    end

end