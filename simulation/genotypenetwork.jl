# Generate the edgelist of the genotype network

using GenoNet.PathwayFramework
using CSV
using Logging

const gns = Genes([1, 2, 3, 4])
const prtns = Proteins([1, 2, 3, 4, 5], 1, 1)
const env = BinaryEnv(prtns, [1], Int[], [5])

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
open("../data/genotype_network.csv", "w") do fp
    edgelist |> CSV.write(fp)
end
