module Utils

using Arpack

export reachable, leadingeigenvec

"""
    _traverse(adjlist::Dict{T, Vector{T}})::Function where {T}
    _traverse(adjlist::Vector{<:Vector{<:Integer}})::Function

Return a void function which traverses a graph via a depth-first search.

The returned void function `trav!(visited, src)` runs a depth-first search on a graph
(represented by its adjacency list `adjlist`) from node `src` and modifies `visited`,
which stores the state of being visited for each node. This function is only for internal
use. See also [`reachable`](@ref).

# Pre-condiction
`visited` must have the same keys/indices as `adjlist`.
"""
function _traverse(adjlist::Dict{T, Vector{T}}) where {T}
    function trav!(visited::Dict{T, Bool}, src::T)::Nothing
        let adjlist = adjlist
            @assert haskey(adjlist, src) "node $src not in the graph"
            visited[src] && return nothing
            visited[src] = true
            foreach(nei -> trav!(visited, nei), adjlist[src])
        end
    end
    return trav!
end

function _traverse(adjlist::Vector{<:Vector{<:Integer}})
    function trav!(visited::BitVector, src::Integer)::Nothing
        let adjlist = adjlist
            @assert checkbounds(Bool, adjlist, src) "node $src not in the graph"
            @inbounds begin
                visited[src] && return nothing
                visited[src] = true
                foreach(nei -> trav!(visited, nei), adjlist[src])
            end
        end
    end
    return trav!
end

"""
    reachable(adjlist::Dict{T, Vector{T}}, srcs::Vector{T})::Dict{T, Bool} where {T}
    reachable(adjlist::Vector{<:Vector{<:Integer}}, srcs::Vector{<:Integer})::BitVector

Return whether nodes in a graph is reachable from a group of source nodes.

For each node in a graph (represented by its adjacency list `adjlist`), determine whether
it is reachable from any of the nodes `srcs`.

# Note
Keys/indices of the adjacency list `adjlist` should contain all the nodes in the graph.

# Examples
```jldoctest
julia> adjlist = Dict(1 => [3,], 2 => [3,], 3 => Int[]);
julia> reachable(adjlist, [2,])
Dict{Int64,Bool} with 3 entries:
  2 => true
  3 => true
  1 => false
```
```jldoctest
julia> adjlist = [[3,], [3,], Int[]];
julia> reachable(adjlist, [2,])
3-element BitArray{1}:
 0
 1
 1
```
"""
function reachable(adjlist::Dict{T, Vector{T}}, srcs::Vector{T}) where {T}
    visited = Dict(node => false for node in keys(adjlist))
    trav! = _traverse(adjlist)
    for src in srcs trav!(visited, src) end
    return visited
end

function reachable(adjlist::Vector{<:Vector{<:Integer}}, srcs::Vector{<:Integer})
    visited = falses(length(adjlist))
    trav! = _traverse(adjlist)
    for src in srcs trav!(visited, src) end
    return visited
end

"""
    leadingeigenvec(mat::AbstractMatrix)::Vector

Return the leading eigenvector of matrix `mat`.

The leading eigenvector of a matrix is the one corresponding to the eigenvalue with the
largest magnitude.
"""
leadingeigenvec(mat::AbstractMatrix) = vec(eigs(mat, nev=1, which=:LM)[2])

end # module Utils
