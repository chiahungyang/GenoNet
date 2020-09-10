module PopGenModel

using ..PathwayFramework
using ..PathwayFramework: _index_activators, _index_products, _abstraction, _default
using Random

export ReproductionMode, reproduce, reproduce!
export ViabilityMode, isviable
export MutationMode, mutate, mutate!
export PopDynamicsMode, evolve
export model, PreAllocPop
export BinaryViability, IdenticalReproductivity, IndependentMutation, ConstantPopSize
export evolve!


# ----------------------------------------------------------------
# Phenotypic responses to the environment

"""
    ViabilityMode

Model assumption on the viability of an individual.

# Implementation
Any custom subtype of `ViabilityMode` should implement an associated method of the
[`isviable`](@ref) function.
"""
abstract type ViabilityMode end

"""
    isviable(pht::AbstractPhenotype, env::AbstractEnv, vm::ViabilityMode)::Bool
    isviable(gt::AbstractGenotype, env::AbstractEnv, vm::ViabilityMode)::Bool

Return whether an individual is viable and survives under the environment.

The individual can be represented by its phenotype `pht` or genotype `gt`, and its
viability under environment `env` is determined following the model specification `vm`.
"""
function isviable end

function isviable(gt::AbstractGenotype, env::AbstractEnv, vm::ViabilityMode)
    @assert iscompatible(gt, env) "inconsistent proteins in genotype and environment"
    return isviable(phenotype(gt, env), env, vm)
end

"""
    BinaryViability <: ViabilityMode

Model specification where viability is definite, and an individual either survives
the selection or not.
"""
struct BinaryViability <: ViabilityMode end

"""
    isviable(pht::BinaryPhenotype{P}, env::BinaryEnv{P}, ::BinaryViability)::Bool

Return whether, in environment `env`, all the essential proteins are present and all the
fatal proteins are absent in phenotype `pht`.
"""
function isviable(pht::BinaryPhenotype{P}, env::BinaryEnv{P}, ::BinaryViability) where {P}
    return all(state(pht, essential(env))) && ! any(state(pht, fatal(env)))
end


# ----------------------------------------------------------------
# Mutation of genotypes

"""
    MutationMode

Model assumption on mutation of the genotype of an individual.

# Implementation
Any custom subtype of `MutationMode` should implement an associated method of the
[`mutate`](@ref) function.
"""
abstract type MutationMode end

"""
    mutate(gt::GT, mm::MutationMode)::GT where {GT <: AbstractGenotype}

Return the individual's genotype after a potential random mutation process from genotype
`gt` under the model assumption `mm`.

Note that the genotype of the individual may remain unchanged.
"""
function mutate end

"""
    mutate!(dst::GT, gt::GT, mm::MutationMode)::Nothing where {GT <: AbstractGenotype}

Modify `dst` to the genotype after a potential random mutation process from genotype `gt`
under the model assumption `mm`.

See also [`mutate`](@ref).
"""
function mutate! end

"""
    IndependentMutation <: MutationMode

Model specification where mutation at each locus is an independent random process.

Mutation at each locus occurs with a constant probability, and it then changes the allele
to a uniformly, randomly chosen sample from other possible alleles.
"""
struct IndependentMutation <: MutationMode
    "per-locus mutation probability"
    prob::Float64
end

"""
    mutate(gt::GT, ::IndependentMutation)::GT where {GT <: DyadicGenotype}

Return the genotype mutated from `gt`, where the expression activators/products of alleles
may be randomly, uniformly sampled from all possible protein activators/products.

See also [`IndependentMutation`](@ref).
"""
function mutate(gt::DyadicGenotype, mm::IndependentMutation)
    actvs, prods = activators(proteins(gt)), products(proteins(gt))
    expr = Dict(g => allele(gt, g) for g in genes(gt))
    for g in Random.randsubseq(genes(gt), mm.prob)
        al = expr[g]
        expr[g] = (Random.rand(actvs) => Random.rand(prods))
        while expr[g] == al global expr[g] = (Random.rand(actvs) => Random.rand(prods)) end
    end
    return DyadicGenotype(genes(gt), proteins(gt), expr)
end

"""
    mutate!(dst::GT, gt::GT, mm::IndependentMutation)::Nothing where {GT <: DyadicGenotype}

Modify `dst` to the genotype mutated from `gt`, where expression activators/products of
alleles may be randomly, uniformly sampled from all possible protein activators/products.
"""
function mutate!(dst::GT, gt::GT, mm::IndependentMutation) where {GT <: DyadicGenotype}
    @assert iscompatible(dst, gt) "inconsistent underlying collection of proteins/genes"
    actvs, prods = _index_activators(proteins(gt)), _index_products(proteins(gt))
    gt_expr, mut_expr = _abstraction(gt), _abstraction(dst)
    for i in eachindex(genes(gt))
        mut_expr[i] = gt_expr[i]
        if Random.rand() < mm.prob
            mut_expr[i] = (Random.rand(actvs) => Random.rand(prods))
            while mut_expr[i] == gt_expr[i]
                global mut_expr[i] = (Random.rand(actvs) => Random.rand(prods))
            end
        end
    end
end


# ----------------------------------------------------------------
# Reproduction and generating offspring from a population

"""
    ReproductionMode

Model assumption on the reproduction process of individuals in a population.

# Implementation
Any custom subtype of `ReproductionMode` should implement an associated method of the
[`reproduce`](@ref) function.
"""
abstract type ReproductionMode end

"""
    reproduce(popl::Vector{GT}, n::Integer, rm::ReproductionMode)::Vector{GT}
        where {GT <: AbstractGenotype}

Return `n` offspring individuals reproduced from population `popl` under the model
assumption `rm`.

Individuals are represented by their genotypes.
"""
function reproduce end

"""
    reproduce!(offsp::Vector{GT}, popl::Vector{GT}, rm::ReproductionMode)::Nothing
        where {GT <: AbstractGenotype}

Populate the offspring population `offsp` with individuals reproduced from population
`popl` under the model assumption `rm`.

Individuals are represented by their genotypes. See also [`reproduce`](@ref).
"""
function reproduce! end

"""
    IdenticalReproductivity <: ReproductionMode

Model specification where individuals in a population all have the same reproductivity.

Namely, offpsring is equally likely to be reproduced by any individual in the population.
"""
struct IdenticalReproductivity <: ReproductionMode end

"""
    reproduce(popl::Vector{GT}, n::Integer, ::IdenticalReproductivity)::Vector{GT}
        where {GT <: AbstractGenotype}

Return `n` offspring individuals that are randomly, uniformly sampled from population
`popl`.

Individuals are represented by their genotypes.
"""
function reproduce(popl::Vector{<:AbstractGenotype}, n::Integer, ::IdenticalReproductivity)
    return map(copy, Random.rand(popl, n))
end

"""
    reproduce!(offsp::Vector{GT}, popl::AbstractVector{GT},
        ::IdenticalReproductivity)::Nothing where {GT <: AbstractGenotype}

Populate the offspring population `offsp` with individuals that are randomly, uniformly
sampled from the population `popl`.
"""
function reproduce!(
    offsp::Vector{GT},
    popl::AbstractVector{GT},
    ::IdenticalReproductivity
    ) where {GT <: AbstractGenotype}

    for i in eachindex(offsp) copy!(offsp[i], Random.rand(popl)) end
end


# ----------------------------------------------------------------
# Dynamics of the population

"""
    PopDynamicsMode

Model assumption on the dynamics of a population.

# Implementation
Any custom subtype of `PopDynamicsMode` should implement an associated method of the
[`evolve`](@ref) function.
"""
abstract type PopDynamicsMode end

"""
    evolve(popl::Vector{GT}, env::AbstractEnv{P}, pdm::PopDynamicsMode,
        specs::Vararg)::Vector{GT} where {GT <: AbstractGenotype{G, P}}

Return the next generation of population `popl` under environment `env`, the model
assumption `pdm` and other model specifications `specs`.

This function is not exported to avoid ambiguity of function names.
"""
function evolve end

"""
    ConstantPopSize <: PopDynamicsMode

Model specification where the population size remains constant.
"""
struct ConstantPopSize <: PopDynamicsMode end

"""
    evolve(popl::Vector{GT}, env::AbstractEnv{P}, ::ConstantPopSize, vm::ViabilityMode,
        rm::ReproductionMode, mm::MutationMode)::Vector{GT} where {G, P,
        GT <: AbstractGenotype{G, P}}

Return the next generation of a fixed-size population that undergoes selection,
reproduction and mutation.

The next generation of population `popl` is evolved under environment `env` and model
specifications `vm`, `rm`, and `mm` for viability, reproduction, and mutation respectively.

# Pre-condition
The subtype `GT <: AbstractGenotype` must support methods of the [`isviable`](@ref),
[`reproduction`](@ref), and [`mutate`](@ref) functions.
"""
function evolve(
    popl::Vector{<:AbstractGenotype{G, P}},
    env::AbstractEnv{P},
    ::ConstantPopSize,
    vm::ViabilityMode,
    rm::ReproductionMode,
    mm::MutationMode
    ) where {G, P}

    survivors = [gt for gt in popl if isviable(gt, env, vm)]
    offspring = reproduce(survivors, length(popl), rm)
    nextgen = [mutate(gt, mm) for gt in offspring]
    return nextgen
end

"""
    evolve!(popl::Vector{GT<:AbstractGenotype{G, P}}, offsp::Vector, isvb::BitVector,
        env::AbstractEnv{P}, ::ConstantPopSize, vm::ViabilityMode, rm::ReproductionMode,
        mm::MutationMode)::Nothing where {G, P}

Evolve a constant-sized population to the next generation under selection, reproduction,
and mutation.

The next generation of population `popl` is evolved under environment `env` and model
specifications `vm`, `rm`, and `mm` for viability, reproduction, and mutation respectively.

The arguments `offsp` and `isvb` are pre-allocated arrays to simulate the evolution more
efficiently.

# Pre-condition
The subtype `GT <: AbstractGenotype` must support methods of the [`isviable`](@ref),
[`reproduction!`](@ref), and [`mutate!`](@ref) functions.
"""
function evolve!(
    popl::Vector{GT},
    offsp::Vector{GT},
    isvb::BitVector,
    env::AbstractEnv{P},
    ::ConstantPopSize,
    vm::ViabilityMode,
    rm::ReproductionMode,
    mm::MutationMode
    ) where {G, P, GT <: AbstractGenotype{G, P}}

    for (i, gt) in Iterators.enumerate(popl) isvb[i] = isviable(gt, env, vm) end
    survivors = @view popl[isvb]
    reproduce!(offsp, survivors, rm)
    for (mut, gt) in Iterators.zip(popl, offsp) mutate!(mut, gt, mm) end
end


# ----------------------------------------------------------------
# Population genetic models

"""
    model(specs::Vararg)::Function

Return a function that evolves a population under model specifications `specs`.
"""
function model end

"""
    model(pdm::PopDynamicsMode, vm::ViabilityMode, rm::ReproductionMode,
        mm:MutationMode)::Function

Return a function that evolves a population with selection, reproduction, and mutation.
"""
function model(pdm::PopDynamicsMode, vm::ViabilityMode, rm::ReproductionMode, mm::MutationMode)
    function _evolve(popl::Vector{<:AbstractGenotype{G, P}}, env::AbstractEnv{P}) where {G, P}
        let pdm = pdm, vm = vm, rm = rm, mm = mm
            return evolve(popl, env, pdm, vm, rm, mm)
        end
    end
    return _evolve
end

"""
    PreAllocPop

Pre-allocating the offspring population to model the evolution more efficiently.
"""
struct PreAllocPop
    "genotype to be allocated"
    GT::Type{<:AbstractGenotype}
    "underlying collection of genes"
    gs::Genes
    "underlying collection of proteins"
    ps::Proteins
    "population size"
    sz::Int
end

"""
    model(::PreAllocPop, ::ConstantPopSize, vm::ViabilityMode, rm::ReproductionMode,
        mm::MutationMode)::Function

Return a function that efficiently evolve a constant-sized population with selection,
reproduction, and mutation.
"""
function model(pap::PreAllocPop, ::ConstantPopSize, vm::ViabilityMode, rm::ReproductionMode, mm::MutationMode)
    offsp = [_default(pap.GT, pap.gs, pap.ps) for i = 1:pap.sz]
    isvb = similar(BitArray, pap.sz)
    function _evolve!(popl::Vector{<:AbstractGenotype{G, P}}, env::AbstractEnv{P}) where {G, P}
        let offsp = offsp, isvb = isvb, vm = vm, rm = rm, mm = mm
            return evolve!(popl, offsp, isvb, env, ConstantPopSize(), vm, rm, mm)
        end
    end
    return _evolve!
end

# ----------------------------------------------------------------

end # module PopGenModel