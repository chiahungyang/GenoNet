# Obtain the genotype with the largest predicted frequency under different constraints

# condition 1: no self-regulating genes and no genes with spare functionality
#              (spareness refers to being activated by an input protein or producing a
#              output that is not part of the environmental condition)
# condition 2: no redundant genes
# condition 3: all genes are activated
# condition 4: no gene is directly activated by a stimulated protein and produces an
#              essential protein

# group 1: all viable genotypes
# group 2: condition 1
# group 3: condition 1 + 2
# group 4: condition 1 + 3
# group 5: condition 1 + 2 + 3
# group 6: condition 2
# group 7: condition 4
# group 8: condition 2 + 4
# group 9: condition 1 + 4
# group 10: condition 1 + 2 + 4
# group 11: condition 1 + 3 + 4
# group 12: condition 1 + 2 + 3 + 4

# This script takes ONE command-line argument, which is the focal case of underlying
# proteins/genes and the environmental condition

using GenoNet.PathwayFramework
using GenoNet.Utils: reachable
using JLD, CSV, DataFrames
using Logging

@info "Loading parameters......"
length(ARGS) == 1 || error("only one command-line argument is accepted")
const CASE, = ARGS

const params = load("../data/$CASE/params.jld")
const gns = Genes(params["gns"]...)
const prtns = Proteins(params["prtns"]...)
const env = BinaryEnv(prtns, params["env"]...)

const gs = collect(gns)
const spare_in = setdiff(input(prtns), stimulated(env))
const spare_out = setdiff(output(prtns), essential(env), fatal(env))

# Load the viable genotypes and their predicted frequencies
@info "Loading data......"
const predict = CSV.File("../data/$CASE/stationary_distribution_apprx.csv") |> DataFrame
const nviable = size(predict, 1)

# Obtain groups of viable genotypes under differnt constraints
# which are represented by their orders in the data
@info "Constructing groups of genotypes......"
const group = Dict("$i" => Int[] for i = 1:12)
for (i, ind) in Iterators.enumerate(predict.gt)
    gt = genotype(DyadicGenotype, gns, prtns, ind)
    als = allele(gt, gs)
    pht = phenotype(gt, env)

    # Check whether constraints are satisfied for different groups
    cond1 = all(actv != prod && actv ∉ spare_in && prod ∉ spare_out for (actv, prod) in als)
    cond2 = length(unique(als)) == length(als)
    cond3 = all(state(pht, actv) for (actv, prod) in als)
    cond4 = !any(actv ∈ stimulated(env) && prod ∈ essential(env) for (actv, prod) in als)

    # Incrementally populate the groups
    push!(group["1"], i)
    cond1 && push!(group["2"], i)
    (cond1 && cond2) && push!(group["3"], i)
    (cond1 && cond3) && push!(group["4"], i)
    (cond1 && cond2 && cond3) && push!(group["5"], i)
    cond2 && push!(group["6"], i)
    cond4 && push!(group["7"], i)
    (cond2 && cond4) && push!(group["8"], i)
    (cond1 && cond4) && push!(group["9"], i)
    (cond1 && cond2 && cond4) && push!(group["10"], i)
    (cond1 && cond3 && cond4) && push!(group["11"], i)
    (cond1 && cond2 && cond3 && cond4) && push!(group["12"], i)
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
ispath("../data/$CASE/prevalentgenotypes/") || mkdir("../data/$CASE/prevalentgenotypes")
for (gp, gt) in pairs(prevalent)
    als = allele(gt, gs)
    actvs, prods = [actv for (actv, prod) in als], [prod for (actv, prod) in als]
    data = DataFrame(g = gs, actv = actvs, prod = prods)
    data |> CSV.write("../data/$CASE/prevalentgenotypes/$gp.csv")
end
