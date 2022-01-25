# Sample viable genotypes in long-term evolution under selection and mutation

# This script takes TWO command-line arguments, where the first one is the focal case of
# the underlying proteins/genes and the environmental condition, and the second one
# indicates the set of parameters for the population genetic simulation

@everywhere using GenoNet.PathwayFramework
@everywhere using GenoNet.PopGenModel
@everywhere using Random
using Distributed, SharedArrays
using JLD, DataFrames, CSV
using Logging

# Hyper-parameters
@info "Loading parameters......"
length(ARGS) == 2 || error("only two command-line arguments are accepted")
const CASE, SIMLT = ARGS

const params = load("../data/$CASE/params.jld")
const gns = Genes(params["gns"]...)
const prtns = Proteins(params["prtns"]...)
const env = BinaryEnv(prtns, params["env"]...)

const simltparams = load("../data/popgensimlt/params/$SIMLT.jld")
const popsz = simltparams["popsz"]  # population size
const mutprob = simltparams["mutprob"]  # per-locus mutation probability

const nexprs = 10000  # number of experiments
const nlins = 1000  # number of lineages per experiment
# number of generations per lineage
const ngens = load("../data/$CASE/convergencetime/$SIMLT.jld", "lowerbound")

# Define the population genetic model
const GT = DyadicGenotype
const preallocation = PreAllocPop(GT, gns, prtns, popsz)
const dynamicsmode = ConstantPopSize()
const viabilitymode = BinaryViability()
const reproductionmode = IdenticalReproductivity()
const mutationmode = IndependentMutation(mutprob)

const _evolve! = model(preallocation, dynamicsmode, viabilitymode, reproductionmode, mutationmode)

# Obtain samples of viable genotypes in long-term evolution
@info "Simulations start"
const samples = SharedArray{Int, 2}(nlins, nexprs)
@sync @distributed for j = 1:nexprs
    anctr = randompopulation(GT, gns, prtns, popsz)  # randomly generated ancestors
    while sum(isviable(gt, env, viabilitymode) for gt in anctr) == 0
        anctr = randompopulation(GT, gns, prtns, popsz)
    end

    popl = defaultpopulation(GT, gns, prtns, popsz)  # pre-allocated population to evolve
    isvb = BitArray(undef, popsz) # pre-allocated indicators of viability

    # Evolve lineages from the ancestral population
    for i = 1:nlins
        for k in eachindex(popl) copy!(popl[k], anctr[k]) end
        for t = 1:ngens _evolve!(popl, env) end
        # Randomly sample the geontype of a viable individual
        for k in eachindex(popl) isvb[k] = isviable(popl[k], env, viabilitymode) end
        samples[i, j] = index(rand(@view popl[isvb]))
    end
    @info "Experiment $j is done"
end

# Output the samples of viable genotypes
@info "Outputing......"
ispath("../data/$CASE/genotypesamples/") || mkdir("../data/$CASE/genotypesamples")
DataFrame(gt = vec(samples)) |> CSV.write("../data/$CASE/genotypesamples/$SIMLT.csv")
