# Estimate the number of generations for the population genetic simulation to converge

# This script takes TWO command-line arguments, where the first one is the focal case of
# the underlying proteins/genes and the environmental condition, and the second one
# indicates the set of parameters for the population genetic simulation

using JLD, DataFrames, CSV
using GenoNet.Utils: variationdist
using GenoNet.PathwayFramework
using Statistics
using Logging

# Load the underlying proteins/genes and the focal per-locus mutation probability
@info "Loading parameters......"
length(ARGS) == 2 || error("only two command-line arguments are accepted")
const CASE, SIMLT = ARGS

const params = load("../data/$CASE/params.jld")
const gns = Genes(params["gns"]...)
const prtns = Proteins(params["prtns"]...)

const simltparams = load("../data/popgensimlt/params/$SIMLTFD.jld")
const mutprob = simltparams["mutprob"]

# Load the number of viable genotypes, the predicted stationary distribution, and its
# empirical distributions
@info "Loading data......"
const viable = CSV.File("../data/$CASE/viable_genotypes.csv") |> DataFrame
const nviable = size(viable, 1)
const apprx = CSV.File("../data/$CASE/stationary_distribution_apprx.csv") |> DataFrame
const distr = CSV.File("../data/$CASE/apprxdistrsamples.csv") |> DataFrame

# Compute the variation distance between the predicted stationary distribution and its
# empirical distributions from finite-sized sampling
@info "Computing variation distances......"
vardist = combine(
    groupby(distr, :expr),
    :freq => (fqs -> variationdist(fqs, apprx.freq)) => :dist
)

# Extending from the theory of Markov chains, compute a lower bound of the convergence time
# with respect to the tolerance of the averaged variation distance
@info "Computing lower bound of convergence time......"
tol = mean(vardist.dist)
nals = length(activators(prtns)) * length(products(prtns))
β = nviable * (mutprob / (nals - 1)) ^ length(gns)
t_cnvg = ceil(Int, log(tol) / log(1 - β))

# Output the estimated convergence time
@info "Outputting......"
ispath("../data/$CASE/convergencetime/") || mkdir("../data/$CASE/convergencetime")
save("../data/$CASE/convergencetime/$SIMLT.jld", "lowerbound", t_cnvg)
