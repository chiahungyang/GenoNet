# Generate the model parameters for different case studies and simulations

using JLD

# ----------------------------------------------------------------
# Parameters for different case studies, including the underlying proteins/genes and the
# environmental condition

# Case "fatal":
# 4 genes, 6 proteins (including 2 input and 2 output ones)
# 1 stimulated and 1 fatal protein
gs = [1, 2, 3, 4]
ps, nin, nout = [1, 2, 3, 4, 5, 6], 2, 2
stml, essnt, ftl = [1], Int[], [6]
gnsparams = [gs]
prtnsparams = [ps, nin, nout]
envparams = [stml, essnt, ftl]
params = Dict("gns" => gnsparams, "prtns" => prtnsparams, "env" => envparams)

case = "fatal"
ispath("../data/$case/") || mkdir("../data/$case")
ispath("../figure/$case/") || mkdir("../figure/$case")
save("../data/$case/params.jld", params)

# Case "essential":
# 4 genes, 6 proteins (including 2 input and 2 output ones)
# 1 stimulated and 1 essential protein
gs = [1, 2, 3, 4]
ps, nin, nout = [1, 2, 3, 4, 5, 6], 2, 2
stml, essnt, ftl = [1], [6], Int[]
gnsparams = [gs]
prtnsparams = [ps, nin, nout]
envparams = [stml, essnt, ftl]
params = Dict("gns" => gnsparams, "prtns" => prtnsparams, "env" => envparams)

case = "essential"
ispath("../data/$case/") || mkdir("../data/$case")
ispath("../figure/$case/") || mkdir("../figure/$case")
save("../data/$case/params.jld", params)

# Case "multiple_fatals":
# 4 genes, 6 proteins (including 2 input and 2 output ones)
# 1 stimulated and 2 fatal proteins
gs = [1, 2, 3, 4]
ps, nin, nout = [1, 2, 3, 4, 5, 6], 2, 2
stml, essnt, ftl = [1], Int[], [5, 6]
gnsparams = [gs]
prtnsparams = [ps, nin, nout]
envparams = [stml, essnt, ftl]
params = Dict("gns" => gnsparams, "prtns" => prtnsparams, "env" => envparams)

case = "multiple_fatals"
ispath("../data/$case/") || mkdir("../data/$case")
ispath("../figure/$case/") || mkdir("../figure/$case")
save("../data/$case/params.jld", params)

# Case "multiple_essentials":
# 4 genes, 6 proteins (including 2 input and 2 output ones)
# 1 stimulated and 2 essential proteins
gs = [1, 2, 3, 4]
ps, nin, nout = [1, 2, 3, 4, 5, 6], 2, 2
stml, essnt, ftl = [1], [5, 6], Int[]
gnsparams = [gs]
prtnsparams = [ps, nin, nout]
envparams = [stml, essnt, ftl]
params = Dict("gns" => gnsparams, "prtns" => prtnsparams, "env" => envparams)

case = "multiple_essentials"
ispath("../data/$case/") || mkdir("../data/$case")
ispath("../figure/$case/") || mkdir("../figure/$case")
save("../data/$case/params.jld", params)

# Case "fatal_multiple_stimuli":
# 4 genes, 6 proteins (including 2 input and 2 output ones)
# 2 stimulated and 1 fatal proteins
gs = [1, 2, 3, 4]
ps, nin, nout = [1, 2, 3, 4, 5, 6], 2, 2
stml, essnt, ftl = [1, 2], Int[], [6]
gnsparams = [gs]
prtnsparams = [ps, nin, nout]
envparams = [stml, essnt, ftl]
params = Dict("gns" => gnsparams, "prtns" => prtnsparams, "env" => envparams)

case = "fatal_multiple_stimuli"
ispath("../data/$case/") || mkdir("../data/$case")
ispath("../figure/$case/") || mkdir("../figure/$case")
save("../data/$case/params.jld", params)

# Case "essential_multiple_stimuli":
# 4 genes, 6 proteins (including 2 input and 2 output ones)
# 2 stimulated and 1 essential proteins
gs = [1, 2, 3, 4]
ps, nin, nout = [1, 2, 3, 4, 5, 6], 2, 2
stml, essnt, ftl = [1, 2], [6], Int[]
gnsparams = [gs]
prtnsparams = [ps, nin, nout]
envparams = [stml, essnt, ftl]
params = Dict("gns" => gnsparams, "prtns" => prtnsparams, "env" => envparams)

case = "essential_multiple_stimuli"
ispath("../data/$case/") || mkdir("../data/$case")
ispath("../figure/$case/") || mkdir("../figure/$case")
save("../data/$case/params.jld", params)

# Case "essential_and_fatal":
# 4 genes, 6 proteins (including 2 input and 2 output ones)
# 1 stimulated, 1 essential and 1 fatal protein
gs = [1, 2, 3, 4]
ps, nin, nout = [1, 2, 3, 4, 5, 6], 2, 2
stml, essnt, ftl = [1], [5], [6]
gnsparams = [gs]
prtnsparams = [ps, nin, nout]
envparams = [stml, essnt, ftl]
params = Dict("gns" => gnsparams, "prtns" => prtnsparams, "env" => envparams)

case = "essential_and_fatal"
ispath("../data/$case/") || mkdir("../data/$case")
ispath("../figure/$case/") || mkdir("../figure/$case")
save("../data/$case/params.jld", params)

# Case "equitable":
# 3 genes, 6 proteins(including 2 input and 2 output ones)
# 2 stimulated, 0 essential and 0 fatal protein
gs = [1, 2, 3]
ps, nin, nout = [1, 2, 3, 4, 5, 6], 2, 2
stml, essnt, ftl = [1, 2], Int[], Int[]
gnsparams = [gs]
prtnsparams = [ps, nin, nout]
envparams = [stml, essnt, ftl]
params = Dict("gns" => gnsparams, "prtns" => prtnsparams, "env" => envparams)

case = "equitable"
ispath("../data/$case/") || mkdir("../data/$case")
ispath("../figure/$case/") || mkdir("../figure/$case")
save("../data/$case/params.jld", params)

# ----------------------------------------------------------------
# Parameters for different population genetic simulations

ispath("../data/popgensimlt/params/") || mkpath("../data/popgensimlt/params")

# Simulation "default":
# 16 individuals, 0.1 per-locus mutation probability
popsz = 16
mutprob = 1e-1
params = Dict("popsz" => popsz, "mutprob" => mutprob)

simlt = "default"
save("../data/popgensimlt/params/$simlt.jld", params)
