
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

smoke.dur.temp <- aggregate(
  cbind(c, N.bar) ~ i + r, data=dat, subset=(w==1), FUN=max
) # used for mean smoking duration, among smokers
smoke.dur <- aggregate(cbind(c, N.bar) ~ r, data=smoke.dur.temp, FUN=mean)

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

est.pois <- coef(m.pois)
var.pois <- diag(vcov(m.pois))

nparam <- length(est.pois)
est.pois <- c(
  est.pois[1:(nparam-2)],
  rep(NA, R-(nparam-2)),
  est.pois[nparam - (1:0)]
)
var.pois <- c(
  var.pois[1:(nparam-2)],
  rep(NA, R-(nparam-2)),
  var.pois[nparam - (1:0)]
)

# vcov.pois <- sandwich::vcovHC(m.pois, type="HC0") # robust variance

est.pois.filename <- "mpois_est.csv"
var.pois.filename <- "mpois_var.csv"

write.table(
  cbind(iter=r, t(est.pois)), file=est.pois.filename,
  append=file.exists(est.pois.filename), quote=F, sep=",",
  row.names=F, col.names=!file.exists(est.pois.filename)
)

write.table(
  cbind(iter=r, t(var.pois)), file=var.pois.filename,
  append=file.exists(var.pois.filename), quote=F, sep=",",
  row.names=F, col.names=!file.exists(var.pois.filename)
)

# saveRDS(m.pois, file=paste0("./models/mpois_iter", r, ".rds"))



print(paste0(Sys.time(), ": Starting Cox analysis"))

## Cox analysis
system.time(m.cox <- coxph(
  Surv(u.prev, u, N.bar) ~ v + c,
  data = dat, subset = (state.prev <= 3),
  method = "breslow"
))
# Takes a few seconds

est.cox <- coef(m.cox)
var.cox <- diag(vcov(m.cox))

est.cox.filename <- "mcox_est.csv"
var.cox.filename <- "mcox_var.csv"

write.table(
  cbind(iter=r, t(est.cox)), file=est.cox.filename,
  append=file.exists(est.cox.filename), quote=F, sep=",",
  row.names=F, col.names=!file.exists(est.cox.filename)
)

write.table(
  cbind(iter=r, t(var.cox)), file=var.cox.filename,
  append=file.exists(var.cox.filename), quote=F, sep=",",
  row.names=F, col.names=!file.exists(var.cox.filename)
)

# saveRDS(m.cox, file=paste0("./models/mcox_iter", r, ".rds"))

