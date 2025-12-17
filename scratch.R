
source("source.R")

###### GENERATE SIMULATION SEEDS ######

set.seed(1)

sim.seeds <- list(.Random.seed)
for (iter in 2:nsim) {
  temp <- runif(10000) # arbitrarily update the random seed
  sim.seeds[[iter]] <- .Random.seed
}
saveRDS(sim.seeds, file=paste0("sim_seeds_nsim", nsim, ".rds"))
