# Plot the stationary distribution of viable genotypes to validate the derivation of the
# the exact solution of the population genetic model

# This script takes TWO command-line arguments, where the first one is the focal case of
# the underlying proteins/genes and the environmental condition, and the second one
# indicates the set of parameters for the population genetic simulation

using JLD, CSV, DataFrames
using Plots, StatsPlots, Measures
using Statistics
using Logging

length(ARGS) == 2 || error("only two command-line arguments are accepted")
const CASE, SIMLT = ARGS

# Set plotting backend to PyPlot
pyplot()

# Load the viable genotypes, the exact stationary distribution from derivation,
# and its confidence band accounting for finite-sized sampling
@info "Loading data......"
const viable = CSV.File("../data/$CASE/viable_genotypes.csv") |> DataFrame
const simltparams = load("../data/popgensimlt/params/$SIMLT.jld")
const mutprob = simltparams["mutprob"]  # per-locus mutation probability
exact = CSV.File("../data/$CASE/stationary_distribution_exact.csv") |> DataFrame
exact = exact[exact.mutprob .== mutprob, Not(:mutprob)]
confintv = CSV.File("../data/$CASE/confidence_interval_exact.csv") |> DataFrame
confintv = confintv[confintv.mutprob .== mutprob, Not(:mutprob)]

# Obtain the ordering of viable genotypes by their predicted frequencies/probabilities
@info "Computing the ordering of genotypes......"
const nviable = size(viable, 1)
const _mean = mean(exact.freq)

sort!(exact, :freq)
exact.order = 1:nviable
const ordering = Dict(r.gt => r.order for r in eachrow(exact))

confintv.order = [ordering[gt] for gt in confintv.gt]
sort!(confintv, :order)

# Load the simulated stationary distribution
@info "Loading the simulated stationary distribution......"
const simlt = CSV.File("../data/$CASE/stationary_distribution_simlt/$SIMLT.csv") |> DataFrame
simlt.order = [ordering[gt] for gt in simlt.gt]

# Plot the confidence band of the derived and the simulation result
@info "Plotting......"
plt = plot(
    xaxis=false,
    xlabel="Viable GRN",
    ylabel="Probability",
    label=""
);

@df simlt scatter!(
    plt,
    :order,
    :freq,
    markersize=3,
    markercolor=colorant"#72DFD7",
    markerstrokewidth=0,
    label=""
);

@df confintv plot!(
    plt,
    :order,
    fill(_mean, nviable),
    ribbon=(_mean .- :lwr, :upr .- _mean),
    linewidth=0,
    fillcolor=colorant"#777777",
    fillalpha=0.5,
    label=""
);

# Output figure
@info "Outputing......"
savefig(
    plot(
        plt,
        left_margin=20mm,
        right_margin=20mm,
        size=(600, 400),
        dpi=300
    ),
    "../figure/$CASE/validate_approximation_extra.png"
)
