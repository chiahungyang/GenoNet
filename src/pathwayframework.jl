module PathwayFramework

export Genes, Proteins, index, input, output, activators, products


# ----------------------------------------------------------------
# Concrete implementation of a collection of genes

"""
Immutable, array-like collection of unique loci with fast access to elements' indices.
"""
struct Genes{T} <: AbstractArray{T, 1}
    "genes/loci"
    gs::Vector{T}
    "indices of loci"
    inds::Dict{T, Int}

    "constructor asserting the uniqueness of loci"
    function Genes{T}(gs, inds) where {T}
        allunique(gs) || throw(DomainError(gs, "loci must be all unique"))
        return new(gs, inds)
    end
end
Genes(gs::Vector{T}, inds::Dict{T, Int}) where {T} = Genes{T}(gs, inds)

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
    index(genes::Genes{T}, g::T)::Int where {T}

Return the index of locus `g` in `genes`.
"""
function index(genes::Genes{T}, g::T)::Int where {T}
    return genes.inds[g]
end


# ----------------------------------------------------------------
# Concrete implementation of a collection of proteins

"""
Immutable, array-like collection of unique proteins with fast access to elements' indices.
"""
struct Proteins{T} <: AbstractArray{T, 1}
    "proteins"
    ps::Vector{T}
    "indices of proteins"
    inds::Dict{T, Int}
    "number of input proteins, which can only be driven by externel stimuli"
    nin::Int
    "number of output proteins, which can only be driven by internal regulation"
    nout::Int

    "constructor asserting the uniqueness of proteins"
    function Proteins{T}(ps, inds, nin, nout) where {T}
        allunique(ps) || throw(DomainError(ps, "proteins must be all unique"))
        return new(ps, inds, nin, nout)
    end
end
function Proteins(ps::Vector{T}, inds::Dict{T, Int}, nin::Int, nout::Int) where {T}
    return Proteins{T}(ps, inds, nin, nout)
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
    Proteins(psin::Vector{T}, psout::Vector{T}, pselse::Vector{T}) where {T}

Construct a Proteins object by concatenating the input, remaining, and output proteins.
"""
function Proteins(psin::Vector{T}, psout::Vector{T}, pselse::Vector{T}) where {T}
    return Proteins([psin; pselse; psout], length(psin), length(psout))
end

"""
    Proteins(actvs::Set{T}, prods::Set{T}) where {T}

Construct a Proteins object from the collections of protein activators and products.

The input, remaining, and output proteins are sorted ascendently and then concatenated
in order.
"""
function Proteins(actvs::Set{T}, prods::Set{T}) where {T}
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
    index(prtns::Proteins{T}, p::T)::Int where {T}

Return the index of protein `p` in `prtns`.
"""
function index(prtns::Proteins{T}, p::T)::Int where {T}
    return prtns.inds[p]
end

"""
    input(prtns::Proteins{T})::Vector{T} where {T}

Return the input proteins in the protein collection `prtns`.
"""
function input(prtns::Proteins{T})::Vector{T} where {T}
    return prtns.ps[begin:begin+prtns.nin-1]
end

"""
    output(prtns::Proteins{T})::Vector{T} where {T}

Return the output proteins in the protein collection `prtns`.
"""
function output(prtns::Proteins{T})::Vector{T} where {T}
    return prtns.ps[end-prtns.nout+1:end]
end

"""
    activators(prtns::Proteins{T})::Vector{T} where {T}

Return those that can serve as expression activators in the protein collection `prtns`.
"""
function activators(prtns::Proteins{T})::Vector{T} where {T}
    return prtns.ps[begin:end-prtns.nout]
end

"""
    products(prtns::Proteins{T})::Vector{T} where {T}

Return those that can serve as expression products in the protein collection `prtns`.
"""
function products(prtns::Proteins{T})::Vector{T} where {T}
    return prtns.ps[begin+prtns.nin:end]
end

# ----------------------------------------------------------------

end # module
