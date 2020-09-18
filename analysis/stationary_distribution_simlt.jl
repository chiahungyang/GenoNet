# Compute the long-term distribution of viable genotypes sampled from simulations

using CSV, DataFrames
using StatsBase
using Logging

const mutprobs = [1e-1, 1e-2, 1e-3]

# Load viable genotypes
@info "Loading data......"
const viable = CSV.File("../data/viable_genotypes.csv") |> DataFrame

# Compute the long-term distribution for simulations with different mutation probabilities
const freqs = Dict{Float64, Vector{Float64}}()
for prob in mutprobs
    # Load samples of genotypes
    @info "Loading data for μ = $prob ......"
    samples = CSV.File("../data/genotypesamples_$prob.csv") |> DataFrame

    # Obtain the  distribution
    @info "Computing distribution for μ = $prob ......"
    prop = proportionmap(samples.gt)
    freqs[prob] = [haskey(prop, gt) ? prop[gt] : 0.0 for gt in viable.ind]
end

# Output the long-term distribution for various mutation probabilities
@info "Building output dataframe......"
distr = DataFrame(gt = Int[], freq = Float64[], mutprob = Float64[])
for prob in mutprobs
    global distr = vcat(distr, DataFrame(gt = viable.ind, freq = freqs[prob], mutprob = prob))
end

@info "Outputing......"
distr |> CSV.write("../data/stationary_distribution_simlt.csv")
