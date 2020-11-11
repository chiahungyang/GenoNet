# Sample genotypes from the stationary distribution under the small-mutation-probability
# approximation, and obtain multiple empirical distributions

using CSV, DataFrames
using StatsBase
using Logging

const nexprs = 1000  # number of empirical distributions to compute
const nsmpls = 10_000_000  # number of samples per empirical distribution

# Load the predicted stationary distribution of genotypes
@info "Loading data......"
const apprx = CSV.File("../data/stationary_distribution_apprx.csv") |> DataFrame
const nviable = size(apprx, 1)

# Pre-allocate to store the samples and the empirical distributions, and pre-compute the
# weights for sampling
@info "Pre-computing......"
const distr = DataFrame(gt = apprx.gt)
const samples = Vector{Int}(undef, nsmpls)
const _counts = Vector{Int}(undef, nviable)

const weights = Weights(apprx.freq)
const _range = 1:nviable

# Obtain samples and compute the empirical distributions
for i = 1:nexprs
    @info "Sampling for experiment $i ......"
    sample!(_range, weights, samples)
    fill!(_counts, 0)
    addcounts!(_counts, samples, _range)
    distr["$i"] = _counts ./ nsmpls
end

# Output the empirical distributions in the tidy format
@info "Constructing dataframe for output......"
data = stack(distr, Not(:gt))
rename!(data, :variable => :expr, :value => :freq)

@info "Outputing......"
data |> CSV.write("../data/apprxdistrsamples.csv")
