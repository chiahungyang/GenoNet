# Generate the indices of all viable genotypes

# This script takes ONE command-line argument, which is the focal case of underlying
# proteins/genes and the environmental condition

using GenoNet.PathwayFramework
using GenoNet.PopGenModel: isviable, BinaryViability
using JLD, CSV
using Logging

@info "Loading parameters......"
length(ARGS) == 1 || error("only one command-line argument is accepted")
const CASE, = ARGS

const params = load("../data/$CASE/params.jld")
const gns = Genes(params["gns"]...)
const prtns = Proteins(params["prtns"]...)
const env = BinaryEnv(prtns, params["env"]...)

@info "Start"

# Iterate over all possible genotypes and examine the viability of each of them
viable = NamedTuple{(:ind,), Tuple{Int}}[]
for (i, gt) in Iterators.enumerate(possiblegenotypes(DyadicGenotype, gns, prtns))
    if isviable(gt, env, BinaryViability()) push!(viable, (ind=index(gt),)) end
    if i % 10000 == 0 @info "$i genotypes explored" end
end

# Output the indices of viable genotypes
@info "Outputing......"
viable |> CSV.write("../data/$CASE/viable_genotypes.csv")
