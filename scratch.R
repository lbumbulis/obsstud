
source("source.R")

###### GENERATE SIMULATION SEEDS ######

RNGkind("L'Ecuyer-CMRG")
set.seed(1999)

sim.seeds <- list(.Random.seed)
for (iter in 2:nsim) {
  sim.seeds[[iter]] <- parallel::nextRNGStream(sim.seeds[[iter-1]])
}
saveRDS(sim.seeds, file=paste0("sim_seeds_nsim", nsim, ".rds"))

RNGkind("default")
