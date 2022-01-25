# Obtain the confidence interval/band of the predicted stationary distribution of viable
# genotypes under the small-mutation-probability approximation

# This script takes ONE command-line argument, which is the focal case of underlying
# proteins/genes and the environmental condition

using CSV, DataFrames
using Statistics
using Logging

const α = 0.95  # confidence level
const lwr = (1 - α) / 2
const upr = 1 - lwr

length(ARGS) == 1 || error("only one command-line argument is accepted")
const CASE, = ARGS

# Load empirical distributions of the predicted stationary distribution
@info "Loading data......"
const distr = CSV.File("../data/$CASE/exactdistrsamples.csv") |> DataFrame

# Compute the confidence interval of the empirical frequencies for each viable genotype
@info "Computing confidence interval......"
confintv = combine(
    groupby(distr, [:gt, :mutprob]),
    :freq => (fqs -> quantile(fqs, lwr)) => :lwr,
    :freq => (fqs -> quantile(fqs, upr)) => :upr
)

# Output the confidence interval/band
@info "Outputing......"
confintv |> CSV.write("../data/$CASE/confidence_interval_exact.csv")
