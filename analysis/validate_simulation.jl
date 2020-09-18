# Plot the stationary distribution of viable genotypes and validate the
# small-mutation-probability approximation (TO BE CHANGED)

using CSV, DataFrames
using VegaLite
using Logging

# Load viable genotypes and the stationary distribution predicted by the approximation and
# obtained from simulations
@info "Loading data......"
const viable = CSV.File("../data/viable_genotypes.csv") |> DataFrame
const nviable, = size(viable)

const apprx = CSV.File("../data/stationary_distribution_apprx.csv") |> DataFrame
const simlt = CSV.File("../data/stationary_distribution_simlt.csv") |> DataFrame

# Obtain the ordering of viable genotypes by their predicted frequencies under the
# approximation
@info "Computing genotype ordering......"
sort!(apprx, :freq)
apprx[:order] = 1:nviable
sort!(apprx, :gt)

const order = Dict(r.gt => r.order for r in eachrow(apprx))
transform!(simlt, :gt => ByRow(gt -> order[gt]) => :order)

# Plot the stationary distribution for the approximation and for the simulations with
# different mutation probability
@info "Constructing data frame for plotting......"
distr = DataFrame(gt = apprx.order, freq = apprx.freq, label = "Approximation")
transform!(simlt, :mutprob => ByRow(p -> "μ = $p") => :label)
distr = vcat(select(simlt, :order => :gt, :freq, :label), distr)

@info "Plotting......"
plt = distr |>
@vlplot(
    :circle,
    x={
        :gt,
        title="Viable genotype"
    },
    y={
        :freq,
        title="Probability"
    },
    color={
        "label:n",
        legend={title="Mutation probability"}
    }
)

# Output the figure
@info "Outputing......"
save("../figure/validate_simulation.pdf", plt)
