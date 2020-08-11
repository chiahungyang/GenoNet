module PathwayFramework

export Genes, index


# ----------------------------------------------------------------
# Concrete implementation of a collection of genes

"""
Immutable, array-like collection of unique loci with fast access to the index of any locus.
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
Base.:(==)(genes::Genes, other::Genes) = genes.gs == other.gs && genes.inds == other.inds

"""
    index(genes::Genes{T}, g::T)::Int where {T}

Return the index of locus `g` in `genes`.
"""
function index(genes::Genes{T}, g::T)::Int where {T}
    return genes.inds[g]
end


# ----------------------------------------------------------------
# Concrete implementation of a collection of proteins

# ----------------------------------------------------------------

end # module
