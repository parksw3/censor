source("seir.R")
source("param.R")

seir_sim <- seir.full(N, beta/N, sigma, gamma, I0, seed=101)

save("seir_sim", file="seir_sim.rda")
