
library(ggplot2)

source("./covariates2/source.R")
source("./covariates2/datagen.R")

###############################################################################
# SAMPLE DATASET
###############################################################################
# sim.seeds <- readRDS("sim_seeds_nsim1000.rds")
# n <- nn
# .Random.seed <- sim.seeds[[1]]
# system.time(dat <- generate.data(n)) # takes < 15min
# saveRDS(dat, "./covariates2/test_dat_iter1.rds", compress="bzip2")

###############################################################################
# MODELLING RESULTS
###############################################################################
#######################################
# Correct model
#######################################
est.cox <- read.csv("./covariates2/sim-results/mcox_cause1_est.csv")
names(est.cox) <- c("iter", paste0("x", c(2, paste0("3.", 1:RB), 4)),"v")
est.cox <- est.cox[order(est.cox$iter),]

# Estimate
round(apply(est.cox[,2:ncol(est.cox)], 2, mean), 3)
# EBias
round(apply(est.cox[,2:ncol(est.cox)], 2, mean) - c(beta1[-(1:2)], gamma1), 3)
# ESE
round(apply(est.cox[,2:ncol(est.cox)], 2, sd), 3)

old.par <- par(mfrow=c(2, 3), mar=c(5,4,1,1))

hist(est.cox$x2, main="", xlab=expression(hat(beta)[12]))
abline(v=c(mean(est.cox$x2), beta1[3]), col=c("black","red"), lty=2, lwd=2)

hist(est.cox$x3.1, main="", xlab=expression(hat(beta)[131]))
abline(v=c(mean(est.cox$x3.1), beta1[4]), col=c("black","red"), lty=2, lwd=2)
hist(est.cox$x3.2, main="", xlab=expression(hat(beta)[132]))
abline(v=c(mean(est.cox$x3.2), beta1[5]), col=c("black","red"), lty=2, lwd=2)
hist(est.cox$x3.3, main="", xlab=expression(hat(beta)[133]))
abline(v=c(mean(est.cox$x3.3), beta1[6]), col=c("black","red"), lty=2, lwd=2)

hist(est.cox$x4, main="", xlab=expression(hat(beta)[14]))
abline(v=c(mean(est.cox$x4), beta1[nvar]), col=c("black","red"), lty=2, lwd=2)

hist(est.cox$v, main="", xlab=expression(hat(gamma)[1]))
abline(v=c(mean(est.cox$v), gamma1), col=c("black","red"), lty=2, lwd=2)
# save as estimates.png with width=750, height=400

par(old.par)

var.cox <- read.csv("./covariates2/sim-results/mcox_cause1_var.csv")
names(var.cox) <- c("iter", paste0("x", c(2, paste0("3.", 1:RB), 4)),"v")
var.cox <- var.cox[order(var.cox$iter),]

# ASE
round(apply(var.cox[,2:ncol(var.cox)], 2, function(x) { mean(sqrt(x)) }), 3)

coverage.cox <- merge(est.cox, var.cox, by="iter", suffixes=c(".est", ".var"))

data.frame(
  param = c(paste0("beta1", c(2, paste0("3.", 1:RB), 4)), "gamma1"),
  coverage = round(c(
    mean(get.coverage(coverage.cox$x2.est, sqrt(coverage.cox$x2.var), beta1[3])), # nsim=1000
    mean(get.coverage(coverage.cox$x3.1.est, sqrt(coverage.cox$x3.1.var), beta1[4])),
    mean(get.coverage(coverage.cox$x3.2.est, sqrt(coverage.cox$x3.2.var), beta1[5])),
    mean(get.coverage(coverage.cox$x3.3.est, sqrt(coverage.cox$x3.3.var), beta1[6])),
    mean(get.coverage(coverage.cox$x4.est, sqrt(coverage.cox$x4.var), beta1[nvar])),
    mean(get.coverage(coverage.cox$v.est, sqrt(coverage.cox$v.var), gamma1))
  ), 3)
)

#######################################
## Incorrect model (omit symptoms)
#######################################
covariate.names <- c(paste0("x", c(2, paste0("3.", 1:RB))),"v")

est.cox <- read.csv("./covariates2/sim-results/mcox_cause1_est_incorrect.csv")
names(est.cox) <- c("iter", covariate.names)
est.cox <- est.cox[order(est.cox$iter),]

# Estimate
round(apply(est.cox[,2:ncol(est.cox)], 2, mean), 3)
# EBias
round(apply(est.cox[,2:ncol(est.cox)], 2, mean) - c(beta1[3:(nvar-1)], gamma1), 3)
# ESE
round(apply(est.cox[,2:ncol(est.cox)], 2, sd), 3)

old.par <- par(mfrow=c(2, 3), mar=c(5,4,1,1))

hist(est.cox$x2, main="", xlab=expression(hat(beta)[12]), xlim=c(beta1[3], max(est.cox$x2)))
abline(v=c(mean(est.cox$x2), beta1[3]), col=c("black","red"), lty=2, lwd=2)

# plot.new()
hist(est.cox$x3.1, main="", xlab=expression(hat(beta)[131]), xlim=c(beta1[4], max(est.cox$x3.1)))
abline(v=c(mean(est.cox$x3.1), beta1[4]), col=c("black","red"), lty=2, lwd=2)
hist(est.cox$x3.2, main="", xlab=expression(hat(beta)[132]))
abline(v=c(mean(est.cox$x3.2), beta1[5]), col=c("black","red"), lty=2, lwd=2)
hist(est.cox$x3.3, main="", xlab=expression(hat(beta)[133]))
abline(v=c(mean(est.cox$x3.3), beta1[6]), col=c("black","red"), lty=2, lwd=2)

hist(est.cox$v, main="", xlab=expression(hat(gamma)[1]))
abline(v=c(mean(est.cox$v), gamma1), col=c("black","red"), lty=2, lwd=2)
par(old.par)

var.cox <- read.csv("./covariates2/sim-results/mcox_cause1_var_incorrect.csv")
names(var.cox) <- c("iter", covariate.names)
var.cox <- var.cox[order(var.cox$iter),]

# ASE
round(apply(var.cox[,2:ncol(var.cox)], 2, function(x) { mean(sqrt(x)) }), 3)

coverage.cox <- merge(est.cox, var.cox, by="iter", suffixes=c(".est", ".var"))

data.frame(
  param = c(paste0("beta1", c(2, paste0("3.", 1:RB))), "gamma1"),
  coverage = c(
    mean(get.coverage(coverage.cox$x2.est, sqrt(coverage.cox$x2.var), beta1[3])), # nsim=1000
    mean(get.coverage(coverage.cox$x3.1.est, sqrt(coverage.cox$x3.1.var), beta1[4])),
    mean(get.coverage(coverage.cox$x3.2.est, sqrt(coverage.cox$x3.2.var), beta1[5])),
    mean(get.coverage(coverage.cox$x3.3.est, sqrt(coverage.cox$x3.3.var), beta1[6])),
    mean(get.coverage(coverage.cox$v.est, sqrt(coverage.cox$v.var), gamma1))
  )
)

var.cox.robust <- read.csv("./covariates2/sim-results/mcox_cause1_var_incorrect_robust.csv")
names(var.cox.robust) <- c("iter", covariate.names)
var.cox.robust <- var.cox.robust[order(var.cox.robust$iter),]

# rASE
round(apply(var.cox.robust[,2:ncol(var.cox.robust)], 2, function(x) { mean(sqrt(x)) }), 3)


