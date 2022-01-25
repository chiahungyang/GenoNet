# Plot the stationary distribution of viable genotypes to validate the
# small-mutation-probability approximation and the convergence of the population genetic
# model

# This script takes TWO command-line arguments, where the first one is the focal case of
# the underlying proteins/genes and the environmental condition, and the second one
# indicates the set of parameters for the population genetic simulation

using CSV, DataFrames
using Plots, StatsPlots, Measures
using Statistics
using Logging

length(ARGS) == 2 || error("only two command-line arguments are accepted")
const CASE, SIMLT = ARGS

# Set plotting backend to PyPlot
pyplot()

# Load the viable genotypes, the stationary distribution predicted by the approximation,
# and its confidence band accounting for finite-sized sampling
@info "Loading data......"
const viable = CSV.File("../data/$CASE/viable_genotypes.csv") |> DataFrame
const apprx = CSV.File("../data/$CASE/stationary_distribution_apprx.csv") |> DataFrame
const confintv = CSV.File("../data/$CASE/confidence_interval_apprx.csv") |> DataFrame

# Obtain the ordering of viable genotypes by their predicted frequencies/probabilities
@info "Computing the ordering of genotypes......"
const nviable = size(viable, 1)
const _mean = mean(apprx.freq)

sort!(apprx, :freq)
apprx.order = 1:nviable
const ordering = Dict(r.gt => r.order for r in eachrow(apprx))

confintv.order = [ordering[gt] for gt in confintv.gt]
sort!(confintv, :order)

# ----------------------------------
# Panel (a):
# Plot the approximated stationary distribution, its confidence band, and the exact
# stationary distributions with respect to differnt per-locus mutation probability

@info "Panel (a):"

# Load the exact stationary distributions calculated from the derivation
@info "Loading the exact stationary distributions......"
const exact = CSV.File("../data/$CASE/stationary_distribution_exact.csv") |> DataFrame
exact.order = [ordering[gt] for gt in exact.gt]

const mutprobs = sort(unique(exact.mutprob), rev=true)
const nprob = length(mutprobs)

# Plot the approximated distribution, its confidence band, and the exact solutions
@info "Plotting......"
plt_a = plot(
    xaxis=false,
    xlabel="Viable GRN",
    ylabel="Probability",
    label="",
    palette=palette([colorant"#EA9197", colorant"#B891EA"], nprob),
    legend=:bottomright
);

@df confintv plot!(
    plt_a,
    :order,
    fill(_mean, nviable),
    ribbon=(_mean .- :lwr, :upr .- _mean),
    linewidth=0,
    fillcolor=colorant"#777777",
    fillalpha=0.3,
    label=""
);

for (i, prob) in Iterators.enumerate(mutprobs)
    @df exact[exact.mutprob .== prob, :] scatter!(
        plt_a,
        :order,
        :freq,
        markersize=3,
        markercolor=i,
        markerstrokewidth=0,
        label="μ = $prob"
    );
end

@df apprx scatter!(
    plt_a,
    :order,
    :freq,
    markersize=3,
    markercolor=colorant"#777777",
    markeralpha=0.3,
    markerstrokewidth=0,
    label=""
);

# ----------------------------------
# Panel (b):
# Plot the confidence band of the approximated distribution and the simulated stationary
# distribution obtaining from numerically evolving the population genetic model

@info "Panel (b):"

# Load the simulated stationary distribution
@info "Loading the simulated stationary distribution......"
const simlt = CSV.File("../data/$CASE/stationary_distribution_simlt/$SIMLT.csv") |> DataFrame
simlt.order = [ordering[gt] for gt in simlt.gt]

# Plot the confidence band of the approximation and the simulation result
@info "Plotting......"
plt_b = plot(
    xaxis=false,
    xlabel="Viable GRN",
    ylabel="Probability",
    label=""
);

@df simlt scatter!(
    plt_b,
    :order,
    :freq,
    markersize=3,
    markercolor=colorant"#72DFD7",
    markerstrokewidth=0,
    label=""
);

@df confintv plot!(
    plt_b,
    :order,
    fill(_mean, nviable),
    ribbon=(_mean .- :lwr, :upr .- _mean),
    linewidth=0,
    fillcolor=colorant"#777777",
    fillalpha=0.5,
    label=""
);

# ----------------------------------
# Combine panels to output

@info "Outputing......"
plt = plot(
    plt_a,
    plt_b,
    title=["(a)" "(b)"],
    titlelocation=:left,
    titlefontsize=10,
    left_margin=20mm,
    right_margin=20mm,
    size=(1200, 400),
    dpi=300
)
savefig(plt, "../figure/$CASE/validate_approximation.png")
