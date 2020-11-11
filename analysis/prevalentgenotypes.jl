# Obtain the genotype with the largest predicted frequency under different constraints

# group 1: all viable genotypes
# group 2: no self-loops
# group 3: no self-loops + no parallel edges
# group 4: no self-loops + all reachable from source
# group 5: no self-loops + no parallel edges + all reachable from source

using GenoNet.PathwayFramework
using GenoNet.Utils: reachable
using CSV, DataFrames
using Logging

const gns = Genes([1, 2, 3, 4])
const prtns = Proteins([1, 2, 3, 4, 5], 1, 1)
const env = BinaryEnv(prtns, [1], Int[], [5])

const gs = collect(gns)

# Load the viable genotypes and their predicted frequencies
@info "Loading data......"
const predict = CSV.File("../data/stationary_distribution_apprx.csv") |> DataFrame
const nviable = size(predict, 1)

# Obtain groups of viable genotypes under differnt constraints
# which are represented by their orders in the data
@info "Constructing groups of genotypes......"
const group = Dict("1" => Int[], "2" => Int[], "3" => Int[], "4" => Int[], "5" => Int[])
for (i, ind) in Iterators.enumerate(predict.gt)
    gt = genotype(DyadicGenotype, gns, prtns, ind)
    als = allele(gt, gs)
    pht = phenotype(gt, env)

    # Check whether constraints are satisfied for different groups
    cond1 = all(actv != prod for (actv, prod) in als)
    cond2 = length(unique(als)) == length(als)
    cond3 = all(state(pht, actv) for (actv, prod) in als)

    # Incrementally populate the groups
    push!(group["1"], i)
    cond1 && push!(group["2"], i)
    (cond1 && cond2) && push!(group["3"], i)
    (cond1 && cond3) && push!(group["4"], i)
    (cond1 && cond2 && cond3) && push!(group["5"], i)
end

# Find the genotype with the largest predicted frequency
@info "Finding prevalent genotypes......"
prevalent = Dict()
for (gp, gts) in pairs(group)
    i = argmax(predict.freq[gts])
    ind = predict.gt[gts[i]]
    prevalent["group$gp"] = genotype(DyadicGenotype, gns, prtns, ind)
end

# Output
@info "Outputing......"
for (gp, gt) in pairs(prevalent)
    als = allele(gt, gs)
    actvs, prods = [actv for (actv, prod) in als], [prod for (actv, prod) in als]
    data = DataFrame(g = gs, actv = actvs, prod = prods)
    data |> CSV.write("../data/prevalentgenotypes/$gp.csv")
end
