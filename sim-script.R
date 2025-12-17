
args <- commandArgs(trailingOnly=TRUE)
r <- as.numeric(args[1])

# library(survival)

source("source.R")
sim.seeds <- readRDS("sim_seeds_nsim1000.rds")



print(paste0(Sys.time(), ": Generating the data"))

.Random.seed <- sim.seeds[[r]]

dat <- generate.data() # for n=1000, usually takes about 30min, sometimes as long as 90min



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
  file = paste0("data_features_iter", r, ".RData")
)


