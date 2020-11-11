# Obtain the confidence interval/band of the predicted stationary distribution of viable
# genotypes under the small-mutation-probability approximation

using CSV, DataFrames
using Statistics
using Logging

const α = 0.95  # confidence level
const lwr = (1 - α) / 2
const upr = 1 - lwr

# Load empirical distributions of the predicted stationary distribution
@info "Loading data......"
distr = CSV.File("../data/apprxdistrsamples.csv") |> DataFrame

# Compute the confidence interval of the empirical frequencies for each viable genotype
@info "Computing confidence interval......"
confintv = combine(
    groupby(distr, :gt),
    :freq => (fqs -> quantile(fqs, lwr)) => :lwr,
    :freq => (fqs -> quantile(fqs, upr)) => :upr
)

# Output the confidence interval/band
@info "Outputing......"
confintv |> CSV.write("../data/confidence_interval_apprx.csv")
