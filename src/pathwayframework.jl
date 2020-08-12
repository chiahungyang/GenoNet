module PathwayFramework

export Genes, Proteins, AbstractPhenotype, BinaryPhenotype
export index, input, output, activators, products, proteins, state


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
    function Genes{G}(gs, inds) where {G}
        allunique(gs) || throw(DomainError(gs, "loci must be all unique"))
        return new(gs, inds)
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
Base.size(genes::Genes) = Base.size(genes.gs)
Base.IndexStyle(::Type{<:Genes}) = Base.IndexLinear()
Base.getindex(genes::Genes, i::Integer) = Base.getindex(genes.gs, i)

# Overwrite and disallow the similar method
Base.similar(genes::Genes) = error("similar not defined for $(typeof(genes))")

# Overwrite the == operator
Base.:(==)(genes::Genes, arr::AbstractArray) = false
function Base.:(==)(genes::Genes, other::Genes)
    return all(name -> getfield(genes, name) == getfield(other, name), fieldnames(Genes))
end

"""
    index(genes::Genes{G}, g::G)::Int where {G}

Return the index of locus `g` in `genes`.
"""
function index(genes::Genes{G}, g::G)::Int where {G}
    return genes.inds[g]
end

"""
    index(genes::Genes{G}, gs::Vector{G})::Vector{Int} where {G}

Return the index of loci `gs` in `genes`.
"""
function index(genes::Genes{G}, gs::Vector{G})::Vector{Int} where {G}
    return [genes.inds[g] for g in gs]
end


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
    function Proteins{P}(ps, inds, nin, nout) where {P}
        allunique(ps) || throw(DomainError(ps, "proteins must be all unique"))
        return new(ps, inds, nin, nout)
    end
end
function Proteins(ps::Vector{P}, inds::Dict{P, Int}, nin::Int, nout::Int) where {P}
    return Proteins{P}(ps, inds, nin, nout)
end

"""
    Proteins(ps::Vector, nin::Int, nout::Int)

Construct a Proteins object in the order of the input, remaining, and output proteins.

The indices of proteins are pre-computed.
"""
function Proteins(ps::Vector, nin::Int, nout::Int)
    return Proteins(ps, Dict(p => i for (i, p) in Iterators.enumerate(ps)), nin, nout)
end

"""
    Proteins(psin::Vector{P}, psout::Vector{P}, pselse::Vector{P}) where {P}

Construct a Proteins object by concatenating the input, remaining, and output proteins.
"""
function Proteins(psin::Vector{P}, psout::Vector{P}, pselse::Vector{P}) where {P}
    return Proteins([psin; pselse; psout], length(psin), length(psout))
end

"""
    Proteins(actvs::Set{P}, prods::Set{P}) where {P}

Construct a Proteins object from the collections of protein activators and products.

The input, remaining, and output proteins are sorted ascendently and then concatenated
in order.
"""
function Proteins(actvs::Set{P}, prods::Set{P}) where {P}
    psin = sort(collect(setdiff(actvs, prods)))
    psout = sort(collect(setdiff(prods, actvs)))
    pselse = sort(collect(intersect(actvs, prods)))
    return Proteins(psin, psout, pselse)
end

# Implement the immutable, array-like interface for Proteins
Base.size(prtns::Proteins) = Base.size(prtns.ps)
Base.IndexStyle(::Type{<:Proteins}) = Base.IndexLinear()
Base.getindex(prtns::Proteins, i::Integer) = Base.getindex(prtns.ps, i)

# Overwrite and disallow the similar method
Base.similar(prtns::Proteins) = error("similar not defined for $(typeof(prtns))")

# Overwrite the == operator
Base.:(==)(prtns::Proteins, arr::AbstractArray) = false
function Base.:(==)(prtns::Proteins, other::Proteins)
    return all(name -> getfield(prtns, name) == getfield(other, name), fieldnames(Proteins))
end

"""
    index(prtns::Proteins{P}, p::P)::Int where {P}

Return the index of protein `p` in `prtns`.
"""
function index(prtns::Proteins{P}, p::P)::Int where {P}
    return prtns.inds[p]
end

"""
    index(prtns::Proteins{P}, p::Vector{P})::Vector{Int} where {P}

Return the indices of proteins `ps` in `prtns`.
"""
function index(prtns::Proteins{P}, ps::Vector{P})::Vector{Int} where {P}
    return [prtns.inds[p] for p in ps]
end

"""
    input(prtns::Proteins{P})::Vector{P} where {P}

Return the input proteins in the protein collection `prtns`.
"""
function input(prtns::Proteins{P})::Vector{P} where {P}
    return prtns.ps[begin:begin+prtns.nin-1]
end

"""
    output(prtns::Proteins{P})::Vector{P} where {P}

Return the output proteins in the protein collection `prtns`.
"""
function output(prtns::Proteins{P})::Vector{P} where {P}
    return prtns.ps[end-prtns.nout+1:end]
end

"""
    activators(prtns::Proteins{P})::Vector{P} where {P}

Return those that can serve as expression activators in the protein collection `prtns`.
"""
function activators(prtns::Proteins{P})::Vector{P} where {P}
    return prtns.ps[begin:end-prtns.nout]
end

"""
    products(prtns::Proteins{P})::Vector{P} where {P}

Return those that can serve as expression products in the protein collection `prtns`.
"""
function products(prtns::Proteins{P})::Vector{P} where {P}
    return prtns.ps[begin+prtns.nin:end]
end


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

Return the collection of underlying proteins of phenotype `pht`.
"""
function proteins end

"""
    state(pht::AbstractPhenotype{P, S}, p::P)::S where {P, S}

Return the chemical state of protein `p` in phenotype `pht`.
"""
function state end

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

function proteins(pht::BinaryPhenotype{P})::Proteins{P} where {P}
    return pht.ps
end

"""
    state(pht::BinaryPhenotype{P}, p::P)::Bool

Return whether protein `p` is present or absent in phenotype `pht`.
"""
function state(pht::BinaryPhenotype{P}, p::P)::Bool where {P}
    return pht.st[index(proteins(pht), p)]
end

"""
    state(pht::BinaryPhenotype{P}, ps::Vector{P})::BitVector

Return whether each protein in `ps` is present or absent in phenotype `pht`.
"""
function state(pht::BinaryPhenotype{P}, ps::Vector{P})::BitVector where {P}
    return pht.st[index(proteins(pht), ps)]
end

# ----------------------------------------------------------------

end # module
