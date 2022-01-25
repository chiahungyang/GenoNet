# Compute the long-term distribution of viable genotypes sampled from simulations

# This script takes TWO command-line arguments, where the first one is the focal case of
# the underlying proteins/genes and the environmental condition, and the second one
# indicates the set of parameters that has been used for the population genetic simulation

using CSV, DataFrames
using StatsBase
using Logging

length(ARGS) == 2 || error("only two command-line arguments are accepted")
const CASE, SIMLT = ARGS

# Load all viable genotypes and the sampled genotypes
@info "Loading data......"
const viable = CSV.File("../data/$CASE/viable_genotypes.csv") |> DataFrame
const samples = CSV.File("../data/$CASE/genotypesamples/$SIMLT.csv") |> DataFrame

# Compute the long-term distribution for simulations with different mutation probabilities
@info "Computing distribution of genotypes ......"
prop = proportionmap(samples.gt)
freqs = [haskey(prop, gt) ? prop[gt] : 0.0 for gt in viable.ind]

# Output the long-term distribution for various mutation probabilities
@info "Outputing......"
distr = DataFrame(gt = viable.ind, freq = freqs)
ispath("../data/$CASE/stationary_distribution_simlt/") || mkdir("../data/$CASE/stationary_distribution_simlt")
distr |> CSV.write("../data/$CASE/stationary_distribution_simlt/$SIMLT.csv")
