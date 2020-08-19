using Test
using GenoNet.PathwayFramework
using GenoNet.PopGenModel

@testset "PopGenModel" begin

    @testset "BinaryViability" begin
        gns = Genes([1, 3, 2])
        prtns = Proteins([1, 3, 2, 4], 1, 1)
        env = BinaryEnv(prtns, [1,], Int[], [4,])
        als_viable = Dict(1 => (1 => 3), 2 => (2 => 3), 3 => (2 => 4))
        gt_viable = DyadicGenotype(gns, prtns, als_viable)
        @test isviable(phenotype(gt_viable, env), env, BinaryViability())
        @test isviable(gt_viable, env, BinaryViability())
        als_inviable = Dict(1 => (1 => 3), 2 => (2 => 3), 3 => (3 => 4))
        gt_inviable = DyadicGenotype(gns, prtns, als_inviable)
        @test ! isviable(phenotype(gt_inviable, env), env, BinaryViability())
        @test ! isviable(gt_inviable, env, BinaryViability())
        wrongenv = BinaryEnv(Proteins([1, 2, 3, 4], 1, 1), [1,], Int[], [4,])
        @test_throws AssertionError isviable(gt_viable, wrongenv, BinaryViability())
    end

end