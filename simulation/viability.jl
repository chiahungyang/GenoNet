# Generate the indices of all viable genotypes

using GenoNet.PathwayFramework
using GenoNet.PopGenModel: isviable, BinaryViability
using CSV
using Logging

const gns = Genes([1, 2, 3, 4])
const prtns = Proteins([1, 2, 3, 4, 5], 1, 1)
const env = BinaryEnv(prtns, [1], Int[], [5])

@info "Start"

# Iterate over all possible genotypes and examine the viability of each of them
viable = NamedTuple{(:ind,), Tuple{Int}}[]
for (i, gt) in Iterators.enumerate(possiblegenotypes(DyadicGenotype, gns, prtns))
    if isviable(gt, env, BinaryViability()) push!(viable, (ind=index(gt),)) end
    if i % 10000 == 0 @info "$i genotypes explored" end
end

# Output the indices of viable genotypes
@info "Outputing......"
open("../data/viable_genotypes.csv", "w") do fp
    viable |> CSV.write(fp)
end
