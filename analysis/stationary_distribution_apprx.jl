# Compute the stationary distribution under the small-mutation-probability approximation,
# which is proportional to the eigenvector centrality in the neutral network

# This script takes ONE command-line argument, which is the focal case of underlying
# proteins/genes and the environmental condition

using CSV
using DataFrames
using SparseArrays
using GenoNet.Utils: leadingeigenvec
using Logging

length(ARGS) == 1 || error("only one command-line argument is accepted")
const CASE, = ARGS

# Load the viable genotypes and the edgelist of the genotype network
@info "Loading data......"
const edgelist = CSV.File("../data/$CASE/genotype_network.csv") |> DataFrame
const viable = CSV.File("../data/$CASE/viable_genotypes.csv") |> DataFrame

const nviable = size(viable, 1)
const ind = Dict(gt => i for (i, gt) in Iterators.enumerate(viable.ind))

# Construct the adjacency matrix of the neutral network and compute its leading eigenvector
@info "Computing approximation......"
const adjmat = spzeros(nviable, nviable)
for e in eachrow(edgelist)
    if haskey(ind, e.src) && haskey(ind, e.trgt)
        adjmat[ind[e.src], ind[e.trgt]] = 1
    end
end
const eigvec = leadingeigenvec(adjmat)

# Output the stationary distribution under the small-mutation-probability approximation
freqs = eigvec ./ sum(eigvec)
distr = DataFrame(gt = viable.ind, freq = freqs)

@info "Outputing......"
distr |> CSV.write("../data/$CASE/stationary_distribution_apprx.csv")
