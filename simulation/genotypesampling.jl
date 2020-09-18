# Sample viable genotypes in long-term evolution under selection and mutation

@everywhere using GenoNet.PathwayFramework
@everywhere using GenoNet.PopGenModel
@everywhere using Random
using Distributed, SharedArrays
using DataFrames, CSV
using Logging

# Hyper-parameters
const gns = Genes([1, 2, 3, 4])
const prtns = Proteins([1, 2, 3, 4, 5], 1, 1)
const popsz = 16  # population size

const env = BinaryEnv(prtns, [1], Int[], [5])
const mutprob = 1e-1  # per-locus mutation probability

const nexprs = 10000  # number of experiments
const nlins = 1000  # number of lineages per experiment
const ngens = round(Int, popsz / mutprob)  # number of generations per lineage

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
DataFrame(gt = vec(samples)) |> CSV.write("../data/genotypesamples_$mutprob.csv")
