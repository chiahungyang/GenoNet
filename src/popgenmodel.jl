module PopGenModel

using ..PathwayFramework
using Random

export ReproductionMode, reproduce!, reproduce
export ViabilityMode, isviable
export BinaryViability, IdenticalReproductivity


# ----------------------------------------------------------------
# Reproduction and generating offspring from a population

"""
    ReproductionMode

Model assumption on the reproduction process of individuals in a population.

# Implementation
Any custom subtype of `ReproductionMode` must implement an associated method of the
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
Any custom subtype of `ViabilityMode` must implement an associated method of the
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

end # module PopGenModel