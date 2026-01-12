
library(ggplot2)

source("./with-symptoms/source.R")

###############################################################################
# MODELLING RESULTS
###############################################################################
est.cox <- read.csv("./with-symptoms/sim-results/mcox_cause2_est.csv")
names(est.cox) <- c("iter","x1","x3","v")
est.cox <- est.cox[order(est.cox$iter),]

hist(est.cox$x1, main="", xlab=expression(hat(beta)[21]))
abline(v=c(mean(est.cox$x1), beta2[1]), col=c("black","red"), lty=2, lwd=2)

hist(est.cox$x3, main="", xlab=expression(hat(beta)[23]))
abline(v=c(mean(est.cox$x3), beta2[3]), col=c("black","red"), lty=2, lwd=2)

hist(est.cox$v, main="", xlab=expression(hat(gamma)[2]))
abline(v=c(mean(est.cox$v), gamma2), col=c("black","red"), lty=2, lwd=2)


var.cox <- read.csv("./with-symptoms/sim-results/mcox_cause2_var.csv")
names(var.cox) <- c("iter","x1","x3","v")
var.cox <- var.cox[order(var.cox$iter),]


get.coverage <- Vectorize(function(est, se, target, conf.level=0.95) {
  crit <- qnorm(conf.level + (1-conf.level)/2)
  return(as.numeric(est - crit*se <= target && target <= est + crit*se))
}, vectorize.args=c("est", "se"))

mean(get.coverage(est.cox$x1, sqrt(var.cox$x1), beta2[1])) # nsim=1000
mean(get.coverage(est.cox$x3, sqrt(var.cox$x3), beta2[3]))
mean(get.coverage(est.cox$v, sqrt(var.cox$v), gamma2))

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


