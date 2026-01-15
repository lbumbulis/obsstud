
library(ggplot2)

source("./with-symptoms/source.R")

###############################################################################
# MODELLING RESULTS
###############################################################################
#######################################
# Correct model
#######################################
est.cox <- read.csv("./with-symptoms/cumulative/sim-results/mcox_cause1_est.csv")
names(est.cox) <- c("iter", paste0("x", c(1:3,5,6)),"v")
est.cox <- est.cox[order(est.cox$iter),]

old.par <- par(mfrow=c(2, 3), mar=c(5,4,1,1))
hist(est.cox$x1, main="", xlab=expression(hat(beta)[11]))
abline(v=c(mean(est.cox$x1), beta1[1]), col=c("black","red"), lty=2, lwd=2)

hist(est.cox$x2, main="", xlab=expression(hat(beta)[12]))
abline(v=c(mean(est.cox$x2), beta1[2]), col=c("black","red"), lty=2, lwd=2)

hist(est.cox$x3, main="", xlab=expression(hat(beta)[13]))
abline(v=c(mean(est.cox$x3), beta1[3]), col=c("black","red"), lty=2, lwd=2)

hist(est.cox$x5, main="", xlab=expression(hat(beta)[15]))
abline(v=c(mean(est.cox$x5), beta1[nvar-1]), col=c("black","red"), lty=2, lwd=2)

hist(est.cox$x6, main="", xlab=expression(hat(beta)[16]))
abline(v=c(mean(est.cox$x6), beta1[nvar]), col=c("black","red"), lty=2, lwd=2)

hist(est.cox$v, main="", xlab=expression(hat(gamma)[1]))
abline(v=c(mean(est.cox$v), gamma1), col=c("black","red"), lty=2, lwd=2)
# save as cause1_cox_estimates.png with width=600, height=500
par(old.par)

var.cox <- read.csv("./with-symptoms/cumulative/sim-results/mcox_cause1_var.csv")
names(var.cox) <- c("iter", paste0("x", c(1:3,5,6)),"v")
var.cox <- var.cox[order(var.cox$iter),]

coverage.cox <- merge(est.cox, var.cox, by="iter", suffixes=c(".est", ".var"))

data.frame(
  param = c(paste0("beta1", c(1:3,5,6)), "gamma1"),
  coverage = c(
    mean(get.coverage(coverage.cox$x1.est, sqrt(coverage.cox$x1.var), beta1[1])), # nsim=1000
    mean(get.coverage(coverage.cox$x2.est, sqrt(coverage.cox$x2.var), beta1[2])),
    mean(get.coverage(coverage.cox$x3.est, sqrt(coverage.cox$x3.var), beta1[3])),
    mean(get.coverage(coverage.cox$x5.est, sqrt(coverage.cox$x5.var), beta1[nvar-1])),
    mean(get.coverage(coverage.cox$x6.est, sqrt(coverage.cox$x6.var), beta1[nvar])),
    mean(get.coverage(coverage.cox$v.est, sqrt(coverage.cox$v.var), gamma1))
  )
)

#######################################
## Incorrect model (omit symptoms)
#######################################
est.cox <- read.csv("./with-symptoms/cumulative/sim-results/mcox_cause1_est_incorrect.csv")
names(est.cox) <- c("iter", paste0("x", c(1:3,5)),"v")
est.cox <- est.cox[order(est.cox$iter),]

old.par <- par(mfrow=c(2, 3), mar=c(5,4,1,1))
hist(est.cox$x1, main="", xlab=expression(hat(beta)[11]), xlim=c(min(est.cox$x1), beta1[1]))
abline(v=c(mean(est.cox$x1), beta1[1]), col=c("black","red"), lty=2, lwd=2)

hist(est.cox$x2, main="", xlab=expression(hat(beta)[12]), xlim=c(beta1[2], max(est.cox$x2)))
abline(v=c(mean(est.cox$x2), beta1[2]), col=c("black","red"), lty=2, lwd=2)

hist(est.cox$x3, main="", xlab=expression(hat(beta)[13]))
abline(v=c(mean(est.cox$x3), beta1[3]), col=c("black","red"), lty=2, lwd=2)

hist(est.cox$x5, main="", xlab=expression(hat(beta)[15]))
abline(v=c(mean(est.cox$x5), beta1[nvar-1]), col=c("black","red"), lty=2, lwd=2)

hist(est.cox$v, main="", xlab=expression(hat(gamma)[1]))
abline(v=c(mean(est.cox$v), gamma1), col=c("black","red"), lty=2, lwd=2)
# save as cause1_cox_estimates.png with width=600, height=500
par(old.par)

var.cox <- read.csv("./with-symptoms/cumulative/sim-results/mcox_cause1_var_incorrect.csv")
names(est.cox) <- c("iter", paste0("x", c(1:3,5)),"v")
var.cox <- var.cox[order(var.cox$iter),]

coverage.cox <- merge(est.cox, var.cox, by="iter", suffixes=c(".est", ".var"))

data.frame(
  param = c(paste0("beta1", c(1,3)), "gamma1"),
  coverage = c(
    mean(get.coverage(coverage.cox$x1.est, sqrt(coverage.cox$x1.var), beta1[1])), # nsim=1000
    mean(get.coverage(coverage.cox$x2.est, sqrt(coverage.cox$x2.var), beta1[2])),
    mean(get.coverage(coverage.cox$x3.est, sqrt(coverage.cox$x3.var), beta1[3])),
    mean(get.coverage(coverage.cox$x5.est, sqrt(coverage.cox$x5.var), beta1[nvar-1])),
    mean(get.coverage(coverage.cox$v.est, sqrt(coverage.cox$v.var), gamma1))
  )
)
