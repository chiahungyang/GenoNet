module PopGenModel

using ..PathwayFramework
using Random

export ReproductionMode, reproduce!, reproduce
export ViabilityMode, isviable
export MutationMode, mutate
export BinaryViability, IdenticalReproductivity, IndependentMutation


# ----------------------------------------------------------------
# Reproduction and generating offspring from a population

"""
    ReproductionMode

Model assumption on the reproduction process of individuals in a population.

# Implementation
Any custom subtype of `ReproductionMode` should implement an associated method of the
[`reproduce!`](@ref) function.
"""
abstract type ReproductionMode end

"""
    reproduce!(offsp::Vector{GT}, popl::Vector{GT}, rm::ReproductionMode)::Nothing where {GT <: AbstractGenotype}

Populate the offspring population `offsp` with individuals reproduced from population
`popl` under the model assumption `rm`.

Individuals are represented by their genotypes. See also [`reproduce`](@ref).
"""
function reproduce! end

"""
    reproduce(popl::Vector{GT}, n::Integer, rm::ReproductionMode)::Vector{GT} where {GT <: AbstractGenotype}

Return `n` offspring individuals reproduced from population `popl` under the model
assumption `rm`.

Individuals are represented by their genotypes.
"""
function reproduce(popl::Vector{GT}, n::Integer, rm::ReproductionMode) where {GT <: AbstractGenotype}
    offsp = similar(popl, n)
    return reproduce!(offsp, popl, rm)
end

"""
    IdenticalReproductivity <: ReproductionMode

Model specification where individuals in a population all have the same reproductivity.

Namely, offpsring is equally likely to be reproduced by any individual in the population.
"""
struct IdenticalReproductivity <: ReproductionMode end

"""
    reproduce!(offsp::Vector{GT}, popl::Vector{GT}, ::IdenticalReproductivity)::Nothing where {GT <: AbstractGenotype}

Populate the offpsring population `offsp` with individuals that are randomly, uniformly
sampled from population `popl`.

Individuals are represented by their genotypes.
"""
function reproduce!(offsp::Vector{GT}, popl::Vector{GT}, ::IdenticalReproductivity) where {GT <: AbstractGenotype}
    Random.rand!(offsp, popl)
end


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
    IndependentMutation <: MutationMode

Model specification where mutation at each locus is an independent random process.

Mutation at each locus occurs with a constant probability `prob`, and it then changes the
the allele to a uniformly, randomly chosen sample from other possible alleles.
"""
struct IndependentMutation <: MutationMode
    prob::Float64
end

"""
    mutate(gt::GT, ::IndependentMutation)::GT where {GT <: DyadicGenotype}

Return the genotype mutated from `gt`, where the expression activators/products of alleles
may be uniformly, randomly sampled from all possible protein activators/products.

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

# ----------------------------------------------------------------

end # module PopGenModel