# Visualize the prevalent genotypes and their corresponding GRNs in different environments

using GenoNet.PathwayFramework
using Plots, GraphRecipes
# using Measures
using JLD, CSV, DataFrames
using StatsBase
using Logging

# ----------------------------------------------------------------
# Utility functions for plotting

"""
    plotgrn_pathway!(plt, gt::DyadicGenotype{G, P}, pos::Dict{P, Tuple{T, T}},
        glabel::Dict{G, String}, plabel::Dict{P, String}, env::BinaryEnv{P})
        where {G, P, T<:Real}

Plot on `plt` the pathway version of the gene regulatory network of genotype `gt`.

The nodes (proteins) are positioned by `pos` and labeled by `plabel`, while the edges
(genes) are labeled by `glabel`. The nodes are colored by their phenotypic state under
environment `env`.
"""
function plotgrn_pathway!(
        plt,
        gt::DyadicGenotype{G, P},
        pos::Dict{P, Tuple{T, T}},
        glabel::Dict{G, String},
        plabel::Dict{P, String},
        env::BinaryEnv{P}
    ) where {G, P, T<:Real}

    gns = genes(gt)
    als = allele(gt, collect(gns))

    # Order nodes (proteins) increasingly by their out-degrees
    outdeg = Dict(p => 0 for p in proteins(gt))
    addcounts!(outdeg, [actv for (actv, prod) in als])
    ps = collect(keys(sort(outdeg, by=k->outdeg[k])))

    # Obtain indices for the nodes, their labels and colors
    _index = Dict(p => i for (i, p) in enumerate(ps))
    nodelaebls = [plabel[p] for p in ps]
    pht = phenotype(gt, env)
    nodecolors = [state(pht, p) ? colorant"#f4c095" : colorant"#afc2d5" for p in ps]

    # Obtain the adjacency list, the edge labels and their width
    adjlist = [Int[] for p in ps]
    edgelabel = Dict()
    ec = Dict()  # counter of multi-edges for a source-target pair
    for (g, (actv, prod)) in zip(gns, als)
        src, trgt = _index[actv], _index[prod]
        push!(adjlist[src], trgt)

        pr = (src, trgt)
        if haskey(ec, pr) ec[pr] += 1 else ec[pr] = 1 end
        edgelabel[(src, trgt, ec[pr])] = glabel[g]
    end

    ew = 2
    edgewidth = (s, d, w) -> ew

    # Handle the weird behavior of GraphRecipes that no node will be plotted when there is
    # no non-self-loop edges
    # Add a non-self-loop edge virtually and hide it with a 0 edge width
    if isempty((src, trgt) for (src, trgts) in enumerate(adjlist) for trgt in trgts if trgt != src)
        @assert length(adjlist) > 1 "current solution can not draw fewer than two nodes"
        push!(adjlist[end], 1)
        edgewidth = (s, d, w) -> (s == d ? ew : 0)
    end

    # Plot the pathway version of the gene regulatory network
    x, y = map(collect, zip([pos[p] for p in ps]...))
    graphplot!(plt, adjlist, x=x, y=y, arrow=arrow(0.6, 0.6), fontsize=20,
               xlims=xlims(plt), ylims=ylims(plt),
               nodeshape=:circle, names=nodelaebls, nodesize=0.4, nodestrokewidth=0,
               nodecolor=nodecolors,
               edgelabel=edgelabel, self_edge_size=0.4, curvature_scalar=0.12,
               edgewidth=edgewidth, edgecolor=[colorant"#777777"], edgelabel_offset=0.15)
end

"""
    plotgrn_conventional!(plt, gt::DyadicGenotype{G, P}, pos::Dict{G, Tuple{T, T}},
        label::Dict{G, String}, env::BinaryEnv{P}) where {G, P, T<:Real}

Plot on `plt` the conventional version of the gene regulatory network of genotype `gt`.

The nodes (genes) are positioned by `pos` and labeled by `label`. The nodes are colored by
their activation status under environment `env`.
"""
function plotgrn_conventional!(
    plt,
    gt::DyadicGenotype{G, P},
    pos::Dict{G, Tuple{T, T}},
    label::Dict{G, String},
    env::BinaryEnv{P}
) where {G, P, T<:Real}

    gns, prtns = genes(gt), proteins(gt)
    als = allele(gt, collect(gns))

    # Obtain the incident edges for each node in the pathway version of the GRN
    incd_out, incd_in = Dict(p => G[] for p in prtns), Dict(p => G[] for p in prtns)
    for (g, (actv, prod)) in zip(gns, als)
        push!(incd_out[actv], g)
        push!(incd_in[prod], g)
    end

    # Order nodes (genes) increasingly by their out-degrees
    outdeg = Dict(g => length(incd_out[prod]) for (g, (actv, prod)) in zip(gns, als))
    gs = collect(keys(sort(outdeg, by=k->outdeg[k])))

    # Obtain indices for the nodes, their labels and colors
    _index = Dict(g => i for (i, g) in enumerate(gs))
    nodelaebls = [label[g] for g in gs]
    pht = phenotype(gt, env)
    activated = Dict(g => state(pht, actv) for (g, (actv, prod)) in zip(gns, als))
    nodecolors = [activated[g] ? colorant"#f4c095" : colorant"#afc2d5" for g in gs]

    # Obtain the adjacency list and the edge widths
    adjlist = [Int[] for g in gs]
    for p in prtns, g1 in incd_in[p], g2 in incd_out[p]
        src, trgt = _index[g1], _index[g2]
        push!(adjlist[src], trgt)
    end

    ew = 2
    edgewidth = (s, d, w) -> ew

    # Handle the weird behavior of GraphRecipes that no node will be plotted when there ss
    # no non-self-loop edges
    # Add a non-self-loop edge virtually and hide it with a 0 edge width
    if isempty((src, trgt) for (src, trgts) in enumerate(adjlist) for trgt in trgts if trgt != src)
        @assert length(adjlist) > 1 "current solution can not draw fewer than two nodes"
        push!(adjlist[end], 1)
        edgewidth = (s, d, w) -> (s == d ? ew : 0)
    end

    # Plot the conventional version of the gene regulatory network
    x, y = map(collect, zip([pos[g] for g in gs]...))
    graphplot!(plt, adjlist, x=x, y=y, arrow=arrow(0.6, 0.6), fontsize=20,
               xlims=xlims(plt), ylims=ylims(plt),
               nodeshape=:rect, names=nodelaebls, nodecolor=nodecolors, nodesize=0.2,
               nodestrokewidth=0,
               self_edge_size=0.6, curvature_scalar=0.12, edgewidth=edgewidth,
               edgecolor=[colorant"#777777"])
end

# ----------------------------------------------------------------
# Main script

const proteinlabel = Dict(1 => "1", 2 => "2", 3 => "3", 4 => "4", 5 => "5", 6 => "6")
const genelabel = Dict(1 => "A", 2 => "B", 3 => "C", 4 => "D")

const rad_p, rad_g, yoffset = 1, 0.6, 1.2
const proteinpos = Dict(
    1 => (rad_p * cos(5 / 6 * 2 * π), rad_p * sin(5 / 6 * 2 * π) + yoffset),
    2 => (rad_p * cos(4 / 6 * 2 * π), rad_p * sin(4 / 6 * 2 * π) + yoffset),
    3 => (rad_p * cos(3 / 6 * 2 * π), rad_p * sin(3 / 6 * 2 * π) + yoffset),
    4 => (rad_p * cos(0 / 6 * 2 * π), rad_p * sin(0 / 6 * 2 * π) + yoffset),
    5 => (rad_p * cos(2 / 6 * 2 * π), rad_p * sin(2 / 6 * 2 * π) + yoffset),
    6 => (rad_p * cos(1 / 6 * 2 * π), rad_p * sin(1 / 6 * 2 * π) + yoffset)
)
const genepos = Dict(
    1 => (rad_g * cos(3 / 4 * 2 * π), rad_g * sin(3 / 4 * 2 * π) - yoffset),
    2 => (rad_g * cos(2 / 4 * 2 * π), rad_g * sin(2 / 4 * 2 * π) - yoffset),
    3 => (rad_g * cos(0 / 4 * 2 * π), rad_g * sin(0 / 4 * 2 * π) - yoffset),
    4 => (rad_g * cos(1 / 4 * 2 * π), rad_g * sin(1 / 4 * 2 * π) - yoffset)
)

# Specify different environments and groups of GRNs under which the prevalent GRN is
# extracted from
const cases_fatal = ["fatal", "multiple_fatals", "fatal_multiple_stimuli"]
const groups_fatal = [1, 2, 3, 4, 5]
const cases_essential = ["essential", "multiple_essentials", "essential_multiple_stimuli"]
const groups_essential = [1, 6, 7, 8]
const cases_both = ["essential_and_fatal"]
const groups_both = [1, 2, 3, 4, 5, 6, 7, 8]

const caselabel = Dict(
    "fatal" => "Environment 1",
    "multiple_fatals" => "Environment 2",
    "fatal_multiple_stimuli" => "Environment 3",
    "essential" => "Environment 4",
    "multiple_essentials" => "Environment 5",
    "essential_multiple_stimuli" => "Environment 6",
    "essential_and_fatal" => "Environment 7"
)
const grouplabel = Dict(
    1 => "Group (\u2170)",
    2 => "Group (\u2171)",
    3 => "Group (\u2172)",
    4 => "Group (\u2173)",
    5 => "Group (\u2174)",
    6 => "Group (\u2175)",
    7 => "Group (\u2176)",
    8 => "Group (\u2177)"
)

# Set plotting backend to PyPlot
pyplot()

# ----------------------------------
# Cases with fatal proteins

plts = []
for (i, case) in enumerate(cases_fatal), (j, group) in enumerate(groups_fatal) 
    # Load the case parameters and the prevalent genotype
    params = load("../data/$case/params.jld")
    gns, prtns = Genes(params["gns"]...), Proteins(params["prtns"]...)
    env = BinaryEnv(prtns, params["env"]...)
    als = CSV.File("../data/$case/prevalentgenotypes/group$group.csv") |> DataFrame
    gt = DyadicGenotype(gns, prtns, Dict(r.g => (r.actv => r.prod) for r in eachrow(als)))

    # Plot the prevalent GRN in both the pathway and the conventional representation
    plt = plot(size=(400, 600), xlims=(-2, 2), ylims=(-3, 3))
    j == 1 && ylabel!(plt, caselabel[case], labelfontsize=36);
    i == 1 && title!(plt, grouplabel[group], titlefontsize=36);
    plotgrn_pathway!(plt, gt, proteinpos, genelabel, proteinlabel, env);
    plotgrn_conventional!(plt, gt, genepos, genelabel, env);
    # plot!(plt, framestyle=:box, bordercolor=colorant"#c4c4c4", margin=0mm);

    push!(plts, plt)
end

@info "Plotting for cases with fatal proteins......"
plt_fatal = plot(
    plts...,
    layout=(length(cases_fatal), length(groups_fatal)),
    size=(400 * length(groups_fatal), 600 * length(cases_fatal)),
    dpi=300
);
@info "Outputting......"
savefig(plt_fatal, "../figure/prevalentgrns_fatal.png")

# ----------------------------------
# Cases with essential proteins

plts = []
for (i, case) in enumerate(cases_essential), (j, group) in enumerate(groups_essential) 
    # Load the case parameters and the prevalent genotype
    params = load("../data/$case/params.jld")
    gns, prtns = Genes(params["gns"]...), Proteins(params["prtns"]...)
    env = BinaryEnv(prtns, params["env"]...)
    als = CSV.File("../data/$case/prevalentgenotypes/group$group.csv") |> DataFrame
    gt = DyadicGenotype(gns, prtns, Dict(r.g => (r.actv => r.prod) for r in eachrow(als)))

    # Plot the prevalent GRN in both the pathway and the conventional representation
    plt = plot(size=(400, 600), xlims=(-2, 2), ylims=(-3, 3))
    j == 1 && ylabel!(plt, caselabel[case], labelfontsize=36);
    i == 1 && title!(plt, grouplabel[group], titlefontsize=36);
    plotgrn_pathway!(plt, gt, proteinpos, genelabel, proteinlabel, env);
    plotgrn_conventional!(plt, gt, genepos, genelabel, env);
    # plot!(plt, framestyle=:box, bordercolor=colorant"#c4c4c4", margin=0mm);

    push!(plts, plt)
end

@info "Plotting for cases with essential proteins......"
plt_essential = plot(
    plts...,
    layout=(length(cases_essential), length(groups_essential)),
    size=(400 * length(groups_essential), 600 * length(cases_essential)),
    dpi=300
);
@info "Outputting......"
savefig(plt_essential, "../figure/prevalentgrns_essential.png")

# ----------------------------------
# Cases with both essential and fatal proteins

plts = []
for (i, case) in enumerate(cases_both), (j, group) in enumerate(groups_both) 
    # Load the case parameters and the prevalent genotype
    params = load("../data/$case/params.jld")
    gns, prtns = Genes(params["gns"]...), Proteins(params["prtns"]...)
    env = BinaryEnv(prtns, params["env"]...)
    als = CSV.File("../data/$case/prevalentgenotypes/group$group.csv") |> DataFrame
    gt = DyadicGenotype(gns, prtns, Dict(r.g => (r.actv => r.prod) for r in eachrow(als)))

    # Plot the prevalent GRN in both the pathway and the conventional representation
    plt = plot(size=(400, 600), xlims=(-2, 2), ylims=(-3, 3))
    j == 1 && ylabel!(plt, caselabel[case], labelfontsize=36);
    i == 1 && title!(plt, grouplabel[group], titlefontsize=36);
    plotgrn_pathway!(plt, gt, proteinpos, genelabel, proteinlabel, env);
    plotgrn_conventional!(plt, gt, genepos, genelabel, env);
    # plot!(plt, framestyle=:box, bordercolor=colorant"#c4c4c4", margin=0mm);

    push!(plts, plt)
end

@info "Plotting for cases with both essential and fatal proteins......"
plt_both = plot(
    plts...,
    layout=(length(cases_both), length(groups_both)),
    size=(400 * length(groups_both), 600 * length(cases_both)),
    dpi=300
);
@info "Outputting......"
savefig(plt_both, "../figure/prevalentgrns_both.png")
