# Sample genotypes from the exact stationary distribution and obtain multiple empirical
# distributions

# This script takes ONE command-line argument, which is the focal case of underlying
# proteins/genes and the environmental condition

using CSV, DataFrames
using StatsBase
using Logging

const nexprs = 1000  # number of empirical distributions to compute
const nsmpls = 10_000_000  # number of samples per empirical distribution

length(ARGS) == 1 || error("only one command-line argument is accepted")
const CASE, = ARGS

# Load the predicted stationary distribution of genotypes
@info "Loading data......"
const exact = CSV.File("../data/$CASE/stationary_distribution_exact.csv") |> DataFrame
# const mutprobs = unique(exact.mutprob)
const mutprobs = [1e-1]
const viable = CSV.File("../data/$CASE/viable_genotypes.csv") |> DataFrame
const nviable = size(viable, 1)

# Pre-allocate to store the samples and the empirical distributions, and pre-compute the
# weights for sampling
@info "Pre-computing......"
const distr = Dict{Float64, DataFrame}()
const samples = Vector{Int}(undef, nsmpls)
const _counts = Vector{Int}(undef, nviable)
const _range = 1:nviable

for prob in mutprobs
    @info "Sampling for μ = $prob ......"
    distr[prob] = DataFrame(gt = exact[exact.mutprob .== prob, :gt])
    weights = Weights(exact[exact.mutprob .== prob, :freq])
    # Obtain samples and compute the empirical distributions
    for i = 1:nexprs
        @info "Sampling for experiment $i ......"
        sample!(_range, weights, samples)
        fill!(_counts, 0)
        addcounts!(_counts, samples, _range)
        distr[prob]["$i"] = _counts ./ nsmpls
    end
end

# Output the empirical distributions in the tidy format
@info "Constructing dataframe for output......"
data = DataFrame(gt = Int[], expr = Int[], freq = Float64[], mutprob = Float64[])
for prob in mutprobs
    stacked = stack(distr[prob], Not(:gt), variable_name=:expr, value_name=:freq)
    stacked.mutprob = prob
    global data = vcat(data, stacked)
end

@info "Outputing......"
data |> CSV.write("../data/$CASE/exactdistrsamples.csv")
