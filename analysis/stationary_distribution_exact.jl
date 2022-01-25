# Compute the stationary distrution derived from the transition matrix with various
# mutation probabilities

# This script takes ONE command-line argument, which is the focal case of underlying
# proteins/genes and the environmental condition

using CSV, DataFrames, JLD
using GenoNet.PathwayFramework
using GenoNet.Utils: leadingeigenvec
using Logging

length(ARGS) == 1 || error("only one command-line argument is accepted")
const CASE, = ARGS

# per-locus mutation probabilities of interest
const mutprobs = [5e-1, 2e-1, 1e-1, 5e-2, 2e-2, 1e-2]

# Load the viable genotypes and the underlying collection of genes/proteins
@info "Loading data......"
const viable = CSV.File("../data/$CASE/viable_genotypes.csv") |> DataFrame
const nviable = size(viable, 1)

const params = load("../data/$CASE/params.jld")
const gns = Genes(params["gns"]...)
const prtns = Proteins(params["prtns"]...)
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
distr |> CSV.write("../data/$CASE/stationary_distribution_exact.csv")
