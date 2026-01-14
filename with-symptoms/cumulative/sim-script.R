
args <- commandArgs(trailingOnly=TRUE)
start.iter <- as.numeric(args[1])
iter.jump <- as.numeric(args[2])
is.correct <- as.logical(args[3])

stop.iter <- start.iter + iter.jump - 1

library(survival)

source("source.R")
source("datagen.R")

sim.seeds <- readRDS("sim_seeds_nsim1000.rds")
.Random.seed <- sim.seeds[[(start.iter %/% iter.jump) + 1]]
# For iter.jump=10: start.iter=1 -> 1st seed, 11 -> 2nd, 21 -> 3rd, ...



for (iter in start.iter:stop.iter) {
  print(paste0(Sys.time(), ": Starting iter ", iter))
  system.time(dat <- generate.data(method="cts", print.increment=1))
  
  if (iter==1) {
    saveRDS(dat, file="test_dat_iter1.rds", compress="bzip2")
  }
  
  # Filter the data to what is needed for analysis
  dat.sub <- dat[-which(dat$state.prev==dat$state & dat$Z %in% 1:2),] # already in absorbing state
  
  print(paste0(Sys.time(), ": Starting Cox analysis"))
  fail.times <- unique(dat.sub$j[which(dat.sub$Z==1)])
  dat.cox <- dat.sub[which(dat.sub$j %in% fail.times),] # only need data from times when someone fails
  
  if (is.correct) {
    m.cox <- coxph(
      Surv(u.prev, u, Z1) ~ factor(E.prev) + x2 + x5 + x6 + v,
      data = dat.cox, method = "breslow"
    )
    coef.est <- coef(m.cox)[c(2,1,3:6)]
    coef.var <- diag(vcov(m.cox))[c(2,1,3:6)]
  } else {
    m.cox <- coxph(
      Surv(u.prev, u, Z1) ~ factor(E.prev) + x2 + x5 + v,
      data = dat.cox, method = "breslow"
    )
    coef.est <- coef(m.cox)[c(2,1,3:5)]
    coef.var <- diag(vcov(m.cox))[c(2,1,3:5)]
  }
  
  coef.est.filename <- paste0("mcox_cause1_est.csv")
  coef.var.filename <- paste0("mcox_cause1_var.csv")
  
  H0.hat <- basehaz(m.cox, centered=F)
  saveRDS(H0.hat, file=paste0("./basehaz/basehaz_cause1_iter", iter, ".rds"), compress="bzip2")
  
  write.table(
    cbind(iter=iter, t(coef.est)), file=coef.est.filename,
    append=file.exists(coef.est.filename), quote=F, sep=",",
    row.names=F, col.names=!file.exists(coef.est.filename)
  )
  
  write.table(
    cbind(iter=iter, t(coef.var)), file=coef.var.filename,
    append=file.exists(coef.var.filename), quote=F, sep=",",
    row.names=F, col.names=!file.exists(coef.var.filename)
  )
}

