module PathwayFramework

using ..Utils: reachable

export Genes, Proteins, index, input, output, activators, products
export AbstractPhenotype, AbstractEnv, AbstractGenotype,
       proteins, state, stimulated, essential, fatal, genes, allele, phenotype, iscompatible
export BinaryPhenotype, BinaryEnv, DyadicGenotype


# ----------------------------------------------------------------
# Concrete implementation of a collection of genes

"""
    Genes{G} <: AbstractArray{G, 1}

Immutable, array-like collection of unique loci of type `G` with fast access to elements'
indices.
"""
struct Genes{G} <: AbstractArray{G, 1}
    "genes/loci"
    gs::Vector{G}
    "indices of loci"
    inds::Dict{G, Int}

    "constructor asserting the uniqueness of loci"
    Genes{G}(gs, inds) where {G} = begin
        allunique(gs) || throw(DomainError(gs, "loci must be all unique"))
        new(gs, inds)
    end
end

Genes(gs::Vector{G}, inds::Dict{G, Int}) where {G} = Genes{G}(gs, inds)

"""
    Genes(gs::Vector)

Construct a Genes object where the indices of loci are pre-computed.
"""
Genes(gs::Vector) = Genes(gs, Dict(g => i for (i, g) in Iterators.enumerate(gs)))

"""
    Genes(gs::Set)

Construct a Genes object where loci are ordered ascendently.
"""
Genes(gs::Set) = Genes(sort(collect(gs)))

# Implement the immutable, array-like interface for Genes
Base.size(gns::Genes) = Base.size(gns.gs)
Base.IndexStyle(::Type{<:Genes}) = Base.IndexLinear()
Base.getindex(gns::Genes, i::Integer) = Base.getindex(gns.gs, i)

# Overwrite and disallow the similar method
Base.similar(::Genes{G}) where {G} = error("similar not defined for Genes{$G}")

# Overwrite the == operator
Base.:(==)(gns::Genes, arr::AbstractArray) = false
function Base.:(==)(gns::Genes{G}, other::Genes{G}) where {G}
    all(name -> getfield(gns, name) == getfield(other, name), fieldnames(Genes))
end

"""
    index(gns::Genes{G}, g::G)::Int where {G}

Return the index of locus `g` in `gns`.
"""
index(gns::Genes{G}, g::G) where {G} = gns.inds[g]

"""
    index(gns::Genes{G}, gs::Vector{G})::Vector{Int} where {G}

Return the indices of loci `gs` in `gns`.
"""
index(gns::Genes{G}, gs::Vector{G}) where {G} = [index(gns, g) for g in gs]


# ----------------------------------------------------------------
# Concrete implementation of a collection of proteins

"""
    Proteins{P} <: AbstractArray{P, 1}

Immutable, array-like collection of unique proteins of type `P` with fast access to
elements' indices.
"""
struct Proteins{P} <: AbstractArray{P, 1}
    "proteins"
    ps::Vector{P}
    "indices of proteins"
    inds::Dict{P, Int}
    "number of input proteins, which can only be driven by external stimuli"
    nin::Int
    "number of output proteins, which can only be driven by internal regulation"
    nout::Int

    "constructor asserting the uniqueness of proteins"
    Proteins{P}(ps, inds, nin, nout) where {P} = begin
        allunique(ps) || throw(DomainError(ps, "proteins must be all unique"))
        new(ps, inds, nin, nout)
    end
end

Proteins(ps::Vector{P}, inds::Dict{P, Int}, nin::Int, nout::Int) where {P} = Proteins{P}(ps, inds, nin, nout)

"""
    Proteins(ps::Vector, nin::Int, nout::Int)

Construct a Proteins object in the order of the input, remaining, and output proteins.

The indices of proteins are pre-computed.
"""
Proteins(ps::Vector, nin::Int, nout::Int) = begin
    nin + nout <= length(ps) || throw(AssertionError("too many input/output proteins"))
    Proteins(ps, Dict(p => i for (i, p) in Iterators.enumerate(ps)), nin, nout)
end

"""
    Proteins(psin::Vector{P}, psout::Vector{P}, pselse::Vector{P}) where {P}

Construct a Proteins object by concatenating the input, remaining, and output proteins.
"""
Proteins(psin::Vector{P}, psout::Vector{P}, pselse::Vector{P}) where {P} = begin
    ps = [psin; pselse; psout]
    Proteins(ps, length(psin), length(psout))
end

"""
    Proteins(actvs::Set{P}, prods::Set{P}) where {P}

Construct a Proteins object from the collections of protein activators and products.

The input, remaining, and output proteins are sorted ascendently and then concatenated
in order.
"""
Proteins(actvs::Set{P}, prods::Set{P}) where {P} = begin
    psin = sort(collect(setdiff(actvs, prods)))
    psout = sort(collect(setdiff(prods, actvs)))
    pselse = sort(collect(intersect(actvs, prods)))
    Proteins(psin, psout, pselse)
end

# Implement the immutable, array-like interface for Proteins
Base.size(prtns::Proteins) = Base.size(prtns.ps)
Base.IndexStyle(::Type{<:Proteins}) = Base.IndexLinear()
Base.getindex(prtns::Proteins, i::Integer) = Base.getindex(prtns.ps, i)

# Overwrite and disallow the similar method
Base.similar(::Proteins{P}) where {P} = error("similar not defined for Proteins{$P}")

# Overwrite the == operator
Base.:(==)(prtns::Proteins, arr::AbstractArray) = false
function Base.:(==)(prtns::Proteins{P}, other::Proteins{P}) where {P}
    all(name -> getfield(prtns, name) == getfield(other, name), fieldnames(Proteins))
end

"""
    index(prtns::Proteins{P}, p::P)::Int where {P}

Return the index of protein `p` in `prtns`.
"""
index(prtns::Proteins{P}, p::P) where {P} = prtns.inds[p]

"""
    index(prtns::Proteins{P}, p::Vector{P})::Vector{Int} where {P}

Return the indices of proteins `ps` in `prtns`.
"""
index(prtns::Proteins{P}, ps::Vector{P}) where {P} = [prtns.inds[p] for p in ps]

"""
    input(prtns::Proteins{P})::Vector{P} where {P}

Return the input proteins in the protein collection `prtns`.
"""
input(prtns::Proteins) = @inbounds prtns.ps[begin:begin+prtns.nin-1]

"""
    output(prtns::Proteins{P})::Vector{P} where {P}

Return the output proteins in the protein collection `prtns`.
"""
output(prtns::Proteins) = @inbounds prtns.ps[end-prtns.nout+1:end]

"""
    activators(prtns::Proteins{P})::Vector{P} where {P}

Return those that can serve as expression activators in the protein collection `prtns`.
"""
activators(prtns::Proteins) = @inbounds prtns.ps[begin:end-prtns.nout]

"""
    products(prtns::Proteins{P})::Vector{P} where {P}

Return those that can serve as expression products in the protein collection `prtns`.
"""
products(prtns::Proteins) = @inbounds prtns.ps[begin+prtns.nin:end]


# ----------------------------------------------------------------
# Abstract interface of a phenotype

"""
    AbstractPhenotype{P, S}

Phenotype as a collective, `S`-typed chemical state of proteins of type `P`.

# Implementation
Any custom subtype of `AbstractPhenotype` should implement methods of the
[`proteins`](@ref) and [`state`](@ref) functions.
"""
abstract type AbstractPhenotype{P, S} end

"""
    proteins(pht::AbstractPhenotype{P, S})::Proteins{P} where {P, S}
    proteins(gt::AbstractGenotype{G, P})::Proteins{P} where {G, P}
    proteins(env::AbstractEnv{P})::Proteins{P} where {P}

Return the collection of underlying proteins.
"""
function proteins end

"""
    state(pht::AbstractPhenotype{P, S}, p::P)::S where {P, S}

Return the chemical state of protein `p` in phenotype `pht`.
"""
function state end

# Implement pretty printing for AbstractPhenotype
Base.show(io::IO, ::MIME"text/plain", pht::AbstractPhenotype) = begin
    print(typeof(pht), ":")
    foreach(p -> print("\n  ", p => state(pht, p)), proteins(pht))
end

# ----------------------------------------------------------------
# Concrete implementation of a phenotype with binary protein states

"""
    BinaryPhenotype{P} <: AbstractPhenotype{P, Bool}

Implementation of a phenotype where the state of a protein of type `P` is binary.

A protein is either present due to external stumili/internal regulation or otherwise
absent.
"""
struct BinaryPhenotype{P} <: AbstractPhenotype{P, Bool}
    "proteins"
    ps::Proteins{P}
    "state of proteins"
    st::BitVector
end

proteins(pht::BinaryPhenotype) = pht.ps

"""
    state(pht::BinaryPhenotype{P}, p::P)::Bool

Return whether protein `p` is present or absent in phenotype `pht`.
"""
state(pht::BinaryPhenotype{P}, p::P) where {P} = @inbounds pht.st[index(proteins(pht), p)]

"""
    state(pht::BinaryPhenotype{P}, ps::Vector{P})::BitVector

Return whether each protein in `ps` is present or absent in phenotype `pht`.
"""
state(pht::BinaryPhenotype{P}, ps::Vector{P}) where {P} = @inbounds pht.st[index(proteins(pht), ps)]


# ----------------------------------------------------------------
# Abstract interface of a selective environment

"""
    AbstractEnv{P}

Selective environment that ecodes information of external stimuli and phenotypic responses
based on proteins of type `P`.

# Implementation
Any custom subtype of `AbstractEnv` should implement a method of the [`proteins`](@ref)
function.
"""
abstract type AbstractEnv{P} end

# ----------------------------------------------------------------
# Conrete implementation of a selective environment acting on the binary state of proteins

"""
    BinaryEnv{P} <: AbstractEnv{P}

Selective environment acting on the binary state of proteins of type `P`.
"""
struct BinaryEnv{P} <: AbstractEnv{P}
    "collection of underlying proteins"
    ps::Proteins{P}
    "stimulated proteins"
    stml::Vector{P}
    "essential proteins"
    essnt::Vector{P}
    "fatal proteins"
    ftl::Vector{P}

    "constructor asserting the consistency of proteins"
    BinaryEnv{P}(ps, stml, essnt, ftl) where {P} = begin
        issubset([stml; essnt; ftl], ps) || throw(AssertionError("proteins not in the collection"))
        isempty(intersect(essnt, ftl)) || throw(AssertionError("no protein is both essential and fatal"))
        issubset(stml, input(ps)) || throw(AssertionError("stimulated protein must be an input protein"))
        issubset([essnt; ftl], output(ps)) || throw(AssertionError("essential/fatal protein must an output protein"))
        new(ps, stml, essnt, ftl)
    end
end

BinaryEnv(ps::Proteins{P}, stml::Vector{P}, essnt::Vector{P}, ftl::Vector{P}) where {P} = BinaryEnv{P}(ps, stml, essnt, ftl)

proteins(env::BinaryEnv) = env.ps

"""
    stimulated(env::BinaryEnv{P})::Vector{P} where {P}

Return proteins of type `P` that are stimulated and thus present in environment `env`.
"""
stimulated(env::BinaryEnv) = env.stml

"""
    essential(env::BinaryEnv{P})::Vector{P} where {P}

Return proteins of type `P` that are required present for viability in environment `env`.
"""
essential(env::BinaryEnv) = env.essnt

"""
    fatal(env::BinaryEnv{P})::Vector{P} where {P}

Return proteins of type `P` that are required absent for viability in environment `env`.
"""
fatal(env::BinaryEnv) = env.ftl


# ----------------------------------------------------------------
# Abstract interface of a genotype

"""
    AbstractGenotype{G, P}

Genotype as collective expression behavior of genes of type `G`, which changes the state of
proteins of type `P`.

# Implementation
Any custom subtype of `AbstractGenotype` should implement methods of the [`genes`](@ref),
[`proteins`](@ref), [`allele`](@ref) and [`phenotype`](@ref) functions.
"""
abstract type AbstractGenotype{G, P} end

"""
    genes(gt::AbstractGenotype{G, P})::Genes{G} where {G, P}

Return the collection of loci of genotype `gt`.
"""
function genes end

"""
    allele(gt::AbstractGenotype{G, P}, g::G) where {G, P}

Return the expression behavior of the allele of locus `g` in genotype `gt`.
"""
function allele end

"""
    phenotype(gt::AbstractGenotype{G, P}, env::AbstractEnv{P})::AbstractPhenotype{P, S}

Return the phenotype of genotype `gt` under selective environment `env`.
"""
function phenotype end

"""
    iscompatible(::AbstractGenotype, ::AbstractEnv)::Bool
    iscompatible(::AbstractGenotype, ::AbstractGenotype)::Bool

Return whether the arguments have the same underlying collection of proteins/genes.
"""
iscompatible(gt::AbstractGenotype, env::AbstractEnv) = begin
    ps_gt, ps_env = proteins(gt), proteins(env)
    ps_gt === ps_env || ps_gt == ps_env
end
iscompatible(gt::AbstractGenotype, other::AbstractGenotype) = begin
    gs_gt, gs_other = genes(gt), genes(other)
    ps_gt, ps_other = proteins(gt), proteins(other)
    (gs_gt === gs_other || gs_gt == gs_other) && (ps_gt === ps_other || ps_gt == ps_other)
end

# ----------------------------------------------------------------
# Concrete implementation of a genotype where expression behavior of a gene is modeled as
# the pair of its presumably sole transcription activator and protein product

"""
    DyadicGenotype{G, P} <: AbstractGenotype{G, P}

Genotype where an allele is modeled as the input-output pair of its expression behavior.

The expression input and output of an allele of a `G`-typed gene is its `P`-typed
transcription activator and protein product respectively, and they are presumably the only
activator-product pair of the allele.
"""
struct DyadicGenotype{G, P} <: AbstractGenotype{G, P}
    "underlying loci"
    gs::Genes{G}
    "underlying proteins"
    ps::Proteins{P}
    "input-output abstraction of expression behavior (denoted via the indices of proteins)"
    abstr::Vector{Pair{Int, Int}}
end

"""
    DyadicGenotype(gs::Genes{G}, ps::Proteins{P}, als::Dict{G, Pair{P, P}})

Construct a genotype from its underly collection of genes `gs` and proteins `ps`, along
with the alleles `als` of `gs` accordingly.
"""
DyadicGenotype(gs::Genes{G}, ps::Proteins{P}, als::Dict{G, Pair{P, P}}) where {G, P} = begin
    # Assert that alleles are consistent with the underlying genes and proteins
    issubset(keys(als), gs) || throw(AssertionError("gene not in the collection"))
    length(als) == length(gs) || throw(AssertionError("alleles of some loci not specified"))
    consistent = all(actv in activators(ps) && prod in products(ps) for (actv, prod) in values(als))
    consistent || throw(AssertionError("protein inconsistent with the collection"))

    DyadicGenotype(gs, ps, [index(ps, als[g].first) => index(ps, als[g].second) for g in gs])
end

genes(gt::DyadicGenotype) = gt.gs
proteins(gt::DyadicGenotype) = gt.ps

"""
    allele(gt::DyadicGenotype{G, P}, g::G)::Pair{P, P}

Return the pair of expression activator and product of the alllele of locus `g`.
"""
allele(gt::DyadicGenotype{G, P}, g::G) where {G, P} = @inbounds begin
    actv, prod = gt.abstr[index(genes(gt), g)]
    proteins(gt)[actv] => proteins(gt)[prod]
end

"""
    allele(gt::DyadicGenotype{G, P}, gs::Vector{G})::Vector{Pair{P, P}}

Return the pair of expression activator and product for each gene in `gs`.
"""
allele(gt::DyadicGenotype{G, P}, gs::Vector{G}) where {G, P} = [allele(gt, g) for g in gs]

"""
    phenotype(gt::DyadicGenotype{G, P}, env::BinaryEnv{P})::BinaryPhenotype{P}

Return the phenotype of genotype `gt` where a protein is present if and only if it falls on
a regulatory pathway triggered by a stimulated protein in environment `env`.
"""
function phenotype(gt::DyadicGenotype{G, P}, env::BinaryEnv{P}) where {G, P}
    @assert iscompatible(gt, env) "inconsistent underlying collection of proteins"
    srcs = [index(proteins(gt), p) for p in stimulated(env)]
    adjlist = [Int[] for p in proteins(gt)]
    for (actv, prod) in gt.abstr append!(adjlist[actv], prod) end
    return BinaryPhenotype(proteins(gt), reachable(adjlist, srcs))
end

# ----------------------------------------------------------------

end # module PathwayFramework
