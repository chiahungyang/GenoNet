module PopGenModel

using ..PathwayFramework

export ViabilityMode, isviable
export BinaryViability


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