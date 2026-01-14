
library(ggplot2)

source("./with-symptoms/source.R")

###############################################################################
# MODELLING RESULTS
###############################################################################
#######################################
# CAUSE 1
#######################################
## Cox
est.cox <- read.csv("./with-symptoms/sim-results/cause1/mcox_cause1_est.csv")
names(est.cox) <- c("iter","x1","x3","x6","v")
est.cox <- est.cox[order(est.cox$iter),]

old.par <- par(mfrow=c(2, 2), mar=c(5,4,1,1))
hist(est.cox$x1, main="", xlab=expression(hat(beta)[11]))
abline(v=c(mean(est.cox$x1), beta1[1]), col=c("black","red"), lty=2, lwd=2)

hist(est.cox$x3, main="", xlab=expression(hat(beta)[13]))
abline(v=c(mean(est.cox$x3), beta1[3]), col=c("black","red"), lty=2, lwd=2)

hist(est.cox$x6, main="", xlab=expression(hat(beta)[16]))
abline(v=c(mean(est.cox$x6), beta1[nvar]), col=c("black","red"), lty=2, lwd=2)

hist(est.cox$v, main="", xlab=expression(hat(gamma)[1]))
abline(v=c(mean(est.cox$v), gamma1), col=c("black","red"), lty=2, lwd=2)
# save as cause1_cox_estimates.png with width=600, height=500
par(old.par)

var.cox <- read.csv("./with-symptoms/sim-results/cause1/mcox_cause1_var.csv")
names(var.cox) <- c("iter","x1","x3","x6","v")
var.cox <- var.cox[order(var.cox$iter),]

coverage.cox <- merge(est.cox, var.cox, by="iter", suffixes=c(".est", ".var"))

data.frame(
  param = c(paste0("beta1", c(1,3,6)), "gamma1"),
  coverage = c(
    mean(get.coverage(coverage.cox$x1.est, sqrt(coverage.cox$x1.var), beta1[1])), # nsim=1000
    mean(get.coverage(coverage.cox$x3.est, sqrt(coverage.cox$x3.var), beta1[3])),
    mean(get.coverage(coverage.cox$x6.est, sqrt(coverage.cox$x6.var), beta1[nvar])),
    mean(get.coverage(coverage.cox$v.est, sqrt(coverage.cox$v.var), gamma1))
  )
)

## Incorrect Cox model (omitting symptom covariate)
est.cox <- read.csv("./with-symptoms/sim-results/cause1/mcox_cause1_est_incorrect.csv")
names(est.cox) <- c("iter","x1","x3","v")
est.cox <- est.cox[order(est.cox$iter),]

old.par <- par(mfrow=c(1, 3), mar=c(5,4,1,1))
hist(est.cox$x1, main="", xlab=expression(hat(beta)[11]))
abline(v=c(mean(est.cox$x1), beta1[1]), col=c("black","red"), lty=2, lwd=2)

hist(est.cox$x3, main="", xlab=expression(hat(beta)[13]))
abline(v=c(mean(est.cox$x3), beta1[3]), col=c("black","red"), lty=2, lwd=2)

hist(est.cox$v, main="", xlab=expression(hat(gamma)[1]))
abline(v=c(mean(est.cox$v), gamma1), col=c("black","red"), lty=2, lwd=2)
# save as cause1_cox_estimates.png with width=600, height=500
par(old.par)

var.cox <- read.csv("./with-symptoms/sim-results/cause1/mcox_cause1_var_incorrect.csv")
names(var.cox) <- c("iter","x1","x3","v")
var.cox <- var.cox[order(var.cox$iter),]

coverage.cox <- merge(est.cox, var.cox, by="iter", suffixes=c(".est", ".var"))

data.frame(
  param = c(paste0("beta1", c(1,3)), "gamma1"),
  coverage = c(
    mean(get.coverage(coverage.cox$x1.est, sqrt(coverage.cox$x1.var), beta1[1])), # nsim=1000
    mean(get.coverage(coverage.cox$x3.est, sqrt(coverage.cox$x3.var), beta1[3])),
    mean(get.coverage(coverage.cox$v.est, sqrt(coverage.cox$v.var), gamma1))
  )
)

#######################################
# CAUSE 2
#######################################
## Cox
est.cox <- read.csv("./with-symptoms/sim-results/mcox_cause2_est_v2.csv")
names(est.cox) <- c("iter","x1","x3","v")
est.cox <- est.cox[order(est.cox$iter),]

par(mfrow=c(1,3))
hist(est.cox$x1, main="", xlab=expression(hat(beta)[21]))
abline(v=c(mean(est.cox$x1), beta2[1]), col=c("black","red"), lty=2, lwd=2)

hist(est.cox$x3, main="", xlab=expression(hat(beta)[23]))
abline(v=c(mean(est.cox$x3), beta2[3]), col=c("black","red"), lty=2, lwd=2)

hist(est.cox$v, main="", xlab=expression(hat(gamma)[2]))
abline(v=c(mean(est.cox$v), gamma2), col=c("black","red"), lty=2, lwd=2)

var.cox <- read.csv("./with-symptoms/sim-results/mcox_cause2_var_v2.csv")
names(var.cox) <- c("iter","x1","x3","v")
var.cox <- var.cox[order(var.cox$iter),]

coverage.cox <- merge(est.cox, var.cox, by="iter", suffixes=c(".est", ".var"))

mean(get.coverage(coverage.cox$x1.est, sqrt(coverage.cox$x1.var), beta2[1])) # nsim=1000
mean(get.coverage(coverage.cox$x3.est, sqrt(coverage.cox$x3.var), beta2[3]))
mean(get.coverage(coverage.cox$v.est, sqrt(coverage.cox$v.var), gamma2))


# Check baseline hazard estimates
basehaz.list <- lapply(1:nsim, function(iter) {
  out <- readRDS(paste0("./with-symptoms/sim-results/raw/basehaz_cause2_iter", iter, ".rds"))
  out$time.idx <- round(out$time * J)
  missing.times <- setdiff(1:((A+tau)*J), out$time.idx)
  missing.low <- missing.times[which(missing.times < min(out$time.idx))]
  missing.high <- missing.times[which(missing.times > max(out$time.idx))]
  missing.mid <- setdiff(missing.times, c(missing.low, missing.high))
  
  out <- rbind(
    out,
    data.frame(hazard=0, time=missing.low/J, time.idx=missing.low),
    data.frame(hazard=NA, time=c(missing.mid/J, missing.high/J), time.idx=c(missing.mid, missing.high))
  )
  out <- out[order(out$time),]
  
  for (tt in missing.mid) {
    tt.idx <- which(out$time.idx==tt)
    out$hazard[tt.idx] <- out$hazard[tt.idx-1]
  }
  return(out$hazard)
})
basehaz <- data.frame(
  hazard = rowMeans(do.call(cbind, basehaz.list)),
  time = (1:((A+tau)*J)) / J
)
plot(hazard ~ time, data=basehaz, type="s")
truehaz <- Vectorize(function(tt) { (rho2*tt)^kap2 })
curve(truehaz, from=0, to=1, col="red", add=T)
# save as basehaz.png with width=500, height=450



## Binomial (logistic)
datagen.method <- "discrete"

est.logis <- read.csv(
  paste0("./with-symptoms/sim-results/", datagen.method, "-datagen/mlogis_cause2_est.csv"),
  fill = T, col.names=c("iter","x1","x3","v", paste0("j", 1:((A+tau)*J)))
)
est.logis <- est.logis[order(est.logis$iter),]

hist(est.logis$x1, main="", xlab=expression(hat(beta)[21]))
abline(v=c(mean(est.logis$x1), beta2[1]), col=c("black","red"), lty=2, lwd=2)

hist(est.logis$x3, main="", xlab=expression(hat(beta)[23]))
abline(v=c(mean(est.logis$x3), beta2[3]), col=c("black","red"), lty=2, lwd=2)

hist(est.logis$v, main="", xlab=expression(hat(gamma)[2]))
abline(v=c(mean(est.logis$v), gamma2), col=c("black","red"), lty=2, lwd=2)

var.logis <- read.csv(
  paste0("./with-symptoms/sim-results/", datagen.method, "-datagen/mlogis_cause2_var.csv"),
  fill = TRUE, col.names = c("iter","x1","x3","v", paste0("j", 1:((A+tau)*J)))
)
var.logis <- var.logis[order(var.logis$iter),]

coverage.logis <- merge(est.logis[,1:4], var.logis[,1:4], by="iter", suffixes=c(".est", ".var"))

mean(get.coverage(coverage.logis$x1.est, sqrt(coverage.logis$x1.var), beta2[1])) # nsim=1000
mean(get.coverage(coverage.logis$x3.est, sqrt(coverage.logis$x3.var), beta2[3]))
mean(get.coverage(coverage.logis$v.est, sqrt(coverage.logis$v.var), gamma2))


