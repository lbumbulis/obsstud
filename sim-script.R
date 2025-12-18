
args <- commandArgs(trailingOnly=TRUE)
r <- as.numeric(args[1])

library(survival)

source("source.R")
sim.seeds <- readRDS("sim_seeds_nsim1000.rds")



print(paste0(Sys.time(), ": Generating the data"))

.Random.seed <- sim.seeds[[r]]

system.time(dat <- generate.data()) # takes < 5min



## Features of the data, to be saved later
v.full <- tapply(dat$v, dat$i, function(x) x[1])
w.full <- tapply(dat$w, dat$i, function(x) x[1])

state.occupancy <- with(dat, table(r, state)) # rows = time points, cols = states; used for plots

transition.count <- with(
  dat, table(state.prev, state, i)
) # used for number of quit attempts (i.e., "2 -> 3" transitions) per person, among smokers

smoke.dur <- aggregate(
  c ~ i, data=dat, FUN=max
) # used for mean smoking duration, among smokers

save(
  v.full, w.full, state.occupancy, transition.count, smoke.dur,
  file = paste0("./data_features/data_features_iter", r, ".RData")
)

if (r==1) {
  saveRDS(dat, file="test_dat_iter1.rds")
}



print(paste0(Sys.time(), ": Starting Poisson analysis"))

## Poisson analysis
system.time(m.pois <- glm(
  N.bar ~ -1 + factor(r) + v + c,
  data = dat,
  subset = (state.prev <= 3),
  family = poisson
))

# est.pois <- coef(m.pois)
# vcov.pois <- sandwich::vcovHC(m.pois, type="HC0") # robust variance

saveRDS(m.pois, file=paste0("./models/mpois_iter", r, ".rds"))



print(paste0(Sys.time(), ": Starting Cox analysis"))

## Cox analysis
system.time(m.cox <- coxph(
  Surv(u.prev, u, N.bar) ~ v + c,
  data = dat, subset = (state.prev <= 3),
  method = "breslow"
))
# Takes a few seconds

# est.cox <- coef(m.cox)
# vcov.cox <- vcov(m.cox)

saveRDS(m.cox, file=paste0("./models/mcox_iter", r, ".rds"))

