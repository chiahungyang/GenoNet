using Test
using GenoNet.PathwayFramework
using GenoNet.PathwayFramework: _default
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
        @test begin
            Random.seed!(12345)
            reproduce(popl, 3, IdenticalReproductivity()) == [gt_inviable, gt_inviable, gt_viable]
        end
        @test begin
            Random.seed!(12345)
            offsp = [_default(DyadicGenotype, gns, prtns) for i = 1:3]
            reproduce!(offsp, popl, IdenticalReproductivity())
            offsp == [gt_inviable, gt_inviable, gt_viable]
        end
    end

    @testset "IndependentMutation" begin
        @test begin
            Random.seed!(12345)
            gt_mut = mutate(gt_inviable, IndependentMutation(0.5))
            gt_mut == DyadicGenotype(gns, prtns, Dict(1 => (1 => 3), 2 => (1 => 2), 3 => (3 => 4)))
        end
        @test begin
            Random.seed!(12345)
            dst = _default(DyadicGenotype, gns, prtns)
            mutate!(dst, gt_inviable, IndependentMutation(0.5))
            dst == DyadicGenotype(gns, prtns, Dict(1 => (1 => 3), 2 => (1 => 2), 3 => (3 => 4)))
        end
        @test_throws AssertionError begin
            Random.seed!(12345)
            dst = _default(DyadicGenotype, Genes([1, 2, 3]), prtns)
            mutate!(dst, gt_inviable, IndependentMutation(0.5))
        end
    end

    @testset "ConstantPopSize" begin
        pdm = ConstantPopSize()
        vm = BinaryViability()
        rm = IdenticalReproductivity()
        mm = IndependentMutation(0.5)
        @test begin
            Random.seed!(12345)
            popl = [gt_viable, gt_inviable, gt_viable]
            nextgen = evolve(popl, env, pdm, vm, rm, mm)
            nextgen == [DyadicGenotype(gns, prtns, Dict(1 => (3 => 2), 2 => (2 => 3), 3 => (3 => 4))),
                        DyadicGenotype(gns, prtns, Dict(1 => (1 => 3), 2 => (2 => 3), 3 => (2 => 4))),
                        DyadicGenotype(gns, prtns, Dict(1 => (1 => 3), 2 => (2 => 3), 3 => (3 => 2)))]
        end
        @test begin
            Random.seed!(12345)
            popl = [DyadicGenotype(gns, prtns, als_viable),
                    DyadicGenotype(gns, prtns, als_inviable),
                    DyadicGenotype(gns, prtns, als_viable)]
            offsp = [_default(DyadicGenotype, gns, prtns) for i = 1:3]
            isvb = similar(BitArray, 3)
            evolve!(popl, offsp, isvb, env, pdm, vm, rm, mm)
            popl == [DyadicGenotype(gns, prtns, Dict(1 => (3 => 2), 2 => (2 => 3), 3 => (3 => 4))),
                     DyadicGenotype(gns, prtns, Dict(1 => (1 => 3), 2 => (2 => 3), 3 => (2 => 4))),
                     DyadicGenotype(gns, prtns, Dict(1 => (1 => 3), 2 => (2 => 3), 3 => (3 => 2)))]
        end
    end

    @testset "Models" begin
        pdm = ConstantPopSize()
        vm = BinaryViability()
        rm = IdenticalReproductivity()
        mm = IndependentMutation(0.5)
        @test begin
            Random.seed!(12345)
            popl = [gt_viable, gt_inviable, gt_viable]
            _evolve = model(pdm, vm, rm, mm)
            nextgen = _evolve(popl, env)
            nextgen == [DyadicGenotype(gns, prtns, Dict(1 => (3 => 2), 2 => (2 => 3), 3 => (3 => 4))),
                        DyadicGenotype(gns, prtns, Dict(1 => (1 => 3), 2 => (2 => 3), 3 => (2 => 4))),
                        DyadicGenotype(gns, prtns, Dict(1 => (1 => 3), 2 => (2 => 3), 3 => (3 => 2)))]
        end
        @test begin
            Random.seed!(12345)
            popl = [DyadicGenotype(gns, prtns, als_viable),
                    DyadicGenotype(gns, prtns, als_inviable),
                    DyadicGenotype(gns, prtns, als_viable)]
            offsp = [_default(DyadicGenotype, gns, prtns) for i = 1:3]
            isvb = similar(BitArray, 3)
            _evolve! = model(PreAllocPop(DyadicGenotype, gns, prtns, 3), pdm, vm, rm, mm)
            _evolve!(popl, env)
            popl == [DyadicGenotype(gns, prtns, Dict(1 => (3 => 2), 2 => (2 => 3), 3 => (3 => 4))),
                     DyadicGenotype(gns, prtns, Dict(1 => (1 => 3), 2 => (2 => 3), 3 => (2 => 4))),
                     DyadicGenotype(gns, prtns, Dict(1 => (1 => 3), 2 => (2 => 3), 3 => (3 => 2)))]
        end
    end

end