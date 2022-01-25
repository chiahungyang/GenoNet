# Generate the edgelist of the genotype network

# This script takes ONE command-line argument, which is the focal case of underlying
# proteins/genes and the environmental condition

using GenoNet.PathwayFramework
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

# Create edgelist by iterate over all genotypes and their mutational neighbors
edgelist = NamedTuple{(:src, :trgt), Tuple{Int, Int}}[]
for (i, gt) in Iterators.enumerate(possiblegenotypes(DyadicGenotype, gns, prtns))
    for mut in possiblemutants(gt)
        push!(edgelist, (src=index(gt), trgt=index(mut)))
    end
    if i % 10000 == 0 @info "$i genotypes explored" end
end

# Output the edgelist of the genotype network
@info "Outputing......"
edgelist |> CSV.write("../data/$CASE/genotype_network.csv")
