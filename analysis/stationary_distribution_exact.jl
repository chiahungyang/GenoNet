# Compute the stationary distrution derived from the transition matrix with various
# mutation probabilities

using CSV, DataFrames
using GenoNet.PathwayFramework
using GenoNet.Utils: leadingeigenvec
using Logging

const mutprobs = [1e-1, 1e-2, 1e-3]  # per-locus probabilities of interest

# Load the viable genotypes and the underlying collection of genes/proteins
@info "Loading data......"
const viable = CSV.File("../data/viable_genotypes.csv") |> DataFrame
const nviable, = size(viable)

const gns = Genes([1, 2, 3, 4])
const prtns = Proteins([1, 2, 3, 4, 5], 1, 1)
const nals = length(activators(prtns)) * length(products(prtns))

# For each per-locus mutation probability of interest, construct the transition matrix and
# compute its leading eigenvector
const eigvec = Dict{Float64, Vector{Float64}}()
const mat = Matrix{Float64}(undef, nviable, nviable)  # pre-allocated transition matrix
for prob in mutprobs
    @info "Constructing transition matrix for μ = $prob ......"
    for (i, ind1) in Iterators.enumerate(viable.ind), (j, ind2) in Iterators.enumerate(viable.ind)
        src = genotype(DyadicGenotype, gns, prtns, ind1)
        dst = genotype(DyadicGenotype, gns, prtns, ind2)
        dist = distance(src, dst)
        mat[i, j] = (prob / (nals - 1)) ^ dist * (1 - prob) ^ (length(gns) - dist)
    end
    @info "Computing solution for μ = $prob ......"
    eigvec[prob] = leadingeigenvec(mat)
end

# Output the stationary distribution for various mutation probability
@info "Building output dataframe......"
distr = DataFrame(gt = Int[], freq = Float64[], mutprob = Float64[])
for prob in mutprobs
    freqs = eigvec[prob] ./ sum(eigvec[prob])
    global distr = vcat(distr, DataFrame(gt = viable.ind, freq = freqs, mutprob = prob))
end

@info "Outputing......"
distr |> CSV.write("../data/stationary_distribution_exact.csv")
