
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
  system.time(dat <- generate.data())
  
  if (iter==1) {
    saveRDS(dat, file="test_dat_iter1.rds", compress="bzip2")
  }
  
  # Filter the data to what is needed for analysis
  dat.sub <- dat[-which(dat$state.prev==dat$state & dat$Z %in% 1:2),] # already in absorbing state
  
  print(paste0(Sys.time(), ": Starting Cox analysis for lung cancer"))
  fail.times <- unique(dat.sub$j[which(dat.sub$Z==1)])
  dat.cox <- dat.sub[which(dat.sub$j %in% fail.times),] # only need data from times when someone fails
  
  if (is.correct) {
    m.cox <- coxph(
      Surv(u.prev, u, Z1) ~ x1 + x2 + x3:B + x4 + v,
      data = dat.cox, method = "breslow"
    )
    coef.order <- c(1:3, (1:RB)+nvar-RB+1, 4:5)
  } else {
    m.cox <- coxph(
      Surv(u.prev, u, Z1) ~ x1 + x2 + x3:B + v,
      data = dat.cox, method = "breslow"
    )
    m.cox.robust <- coxph(
      Surv(u.prev, u, Z1) ~ x1 + x2 + x3:B + v + cluster(i),
      data = dat.cox, method = "breslow", robust=T
    )
    coef.order <- c(1:3, (1:RB)+nvar-RB, 4)
  }
  coef.est <- coef(m.cox)[coef.order]
  coef.var <- diag(vcov(m.cox))[coef.order]
  
  suffix <- ifelse(is.correct, "", "_incorrect")
  coef.est.filename <- paste0("mcox_cause1_est", suffix, ".csv")
  coef.var.filename <- paste0("mcox_cause1_var", suffix, ".csv")
  
  # H0.hat <- basehaz(m.cox, centered=F)
  # saveRDS(H0.hat, file=paste0("./basehaz/basehaz_cause1_iter", iter, ".rds"), compress="bzip2")
  
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
  
  if (!is.correct) {
    coef.var.robust <- diag(vcov(m.cox.robust))[coef.order]
    coef.var.robust.filename <- paste0("mcox_cause1_var_incorrect_robust.csv")
    
    write.table(
      cbind(iter=iter, t(coef.var.robust)), file=coef.var.robust.filename,
      append=file.exists(coef.var.robust.filename), quote=F, sep=",",
      row.names=F, col.names=!file.exists(coef.var.robust.filename)
    )
  }
  
  # if (is.correct) {
  #   print(paste0(Sys.time(), ": Starting Cox analysis for lung cancer symptoms"))
  #   dat.sub <- dat[which(dat$Z.prev==0),]
  #   fail.times <- unique(dat.sub$j[which(dat.sub$Z==10)])
  #   dat.cox <- dat.sub[which(dat.sub$j %in% fail.times),] # only need data from times when someone fails
  #   
  #   m.cox <- coxph(
  #     Surv(u.prev, u, Z10) ~ x1 + x2 + x3:B + v,
  #     data = dat.cox, method = "breslow"
  #   )
  #   coef.order <- c(1:3, (1:RB)+nvar-RB, 4)
  #   coef.est <- coef(m.cox)[coef.order]
  #   coef.var <- diag(vcov(m.cox))[coef.order]
  #   
  #   coef.est.filename <- "mcox_cause1p_est.csv"
  #   coef.var.filename <- "mcox_cause1p_var.csv"
  #   
  #   write.table(
  #     cbind(iter=iter, t(coef.est)), file=coef.est.filename,
  #     append=file.exists(coef.est.filename), quote=F, sep=",",
  #     row.names=F, col.names=!file.exists(coef.est.filename)
  #   )
  #   
  #   write.table(
  #     cbind(iter=iter, t(coef.var)), file=coef.var.filename,
  #     append=file.exists(coef.var.filename), quote=F, sep=",",
  #     row.names=F, col.names=!file.exists(coef.var.filename)
  #   )
  # }
}

