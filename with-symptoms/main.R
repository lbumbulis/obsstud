
library(ggplot2)

source("source.R")

###############################################################################
# DATA FEATURES
###############################################################################
# attempt <- 6 # discrete-time fast data generation
# attempt <- 7 # slow data generation
attempt <- 8 # continuous-time fast data generation

###### State occupancy ################
all.state.occupancy <- lapply(1:nsim, function(iter) {
  filename <- paste0("./sim-results/attempt", attempt, "/data_features_iter", iter, ".RData")
  if (file.exists(filename)) {
    load(filename)
    return(state.occupancy)
  } else {
    print(paste0("Warning: file for iter ", iter, " does not exist"))
    return(NULL)
  }
})
all.state.occupancy <- Filter(Negate(is.null), all.state.occupancy)

plot.dat <- as.data.frame(Reduce("+", all.state.occupancy) / length(all.state.occupancy) / nn)
plot.dat$r <- as.numeric(as.character(plot.dat$r))
plot.dat$is.cancer <- plot.dat$state %in% 4:6

lwd <- 1.5
cols <- c("black","red","blue","green","orange","purple")
ggplot(plot.dat, aes(x=r, y=Freq, group=state, col=state, linetype=is.cancer)) +
  geom_line() +
  ylab("state occupancy prob.") +
  scale_colour_manual(
    labels = c(
      expression("(" * 0^o * ", 0)"), "(1, 0)", "(0, 0)",
      expression("(" * 0^o * ", 1)"), "(1, 1)", "(0, 1)"
    ),
    values = cols
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.key.size=unit(1.5, "lines")
  )
# state_occupancy.png; width=500, height=350

###### Quit attempts ##################
all.transition.count <- lapply(1:nsim, function(iter) {
  filename <- paste0("./sim-results/attempt", attempt, "/data_features_iter", iter, ".RData")
  if (file.exists(filename)) {
    load(filename)
    return(transition.count[,,which(w.full==1)])
  } else {
    print(paste0("Warning: file for iter ", iter, " does not exist"))
    return(NULL)
  }
})
all.transition.count <- Filter(Negate(is.null), all.transition.count)

# Mean number of quit attempts among ever-smokers
quit.dat <- sapply(1:length(all.transition.count), function(iter) {
  mean(sapply(1:(dim(all.transition.count[[iter]])[3]), function(i) {
    mat <- all.transition.count[[iter]][,,i]
    if (mat[1,2]==0) { # never started smoking
      # print(paste0("iter=", iter, ", i=", i, " never started smoking"))
      return(NA)
    } else {
      return(mat[2,3])
    }
  }), na.rm=T)
})
hist(quit.dat, xlab="mean # quit attempts", main="") # nquit_hist.png; width=500, height=400
mean(quit.dat) # attempt 6: 0.922, attempt 7: 0.915, attempt 8: 0.919

###### Smoking duration ###############
all.smoke.dur <- lapply(1:nsim, function(iter) {
  filename <- paste0("./sim-results/attempt", attempt, "/data_features_iter", iter, ".RData")
  if (file.exists(filename)) {
    load(filename)
    return(smoke.dur)
  } else {
    print(paste0("Warning: file for iter ", iter, " does not exist"))
    return(NULL)
  }
})
all.smoke.dur <- Filter(Negate(is.null), all.smoke.dur)

if (attempt <= 7) {
  # Mean (pre-failure) smoking duration among ever-smokers
  dur.dat <- sapply(1:length(all.smoke.dur), function(iter) {
    temp <- all.smoke.dur[[iter]]
    mean(temp$c[which(temp$c > 0)])
  })
  hist(dur.dat, xlab="mean smoking duration", main="") # smokedur_hist.png; width=500, height=400
  mean(dur.dat) # attempt 6: 0.332, attempt 7: 0.330
} else {
  plot.dat2 <- as.data.frame(Reduce("+", all.smoke.dur) / length(all.smoke.dur))
  
  plot(
    c ~ r, data=plot.dat2, type="l", ylim=c(0,1), ylab="",
    main = expression(atop("Empirical mean cumulative exposure", "and failure probability (w = 1)"))
  )
  lines(N.bar ~ r, data=plot.dat2, col="red")
  legend(
    "topleft", lty=1, col=c("black","red"),
    legend = c(expression("Mean" ~ c[r]), "Failure probability")
  )
  
  plot.dat2$c[which(plot.dat2$r==R)] # attempt 8: 0.327
}

###############################################################################
# MODELLING RESULTS
###############################################################################
# attempt <- 6 # discrete-time fast data generation
# attempt <- 7 # slow data generation
attempt <- 8 # continuous-time fast data generation

xlimits <- c(-0.05, 0.45)
ylimits <- c(2, 3.6)

## Poisson
thetahat.pois <- read.csv(paste0("./sim-results/attempt", attempt, "/mpois_est.csv"))
thetahat.pois <- thetahat.pois[, c(1, ncol(thetahat.pois) - 1:0)]
names(thetahat.pois) <- c("iter","alpha1","beta1")

plot(beta1 ~ alpha1, data=thetahat.pois, xlim=xlimits, ylim=ylimits)
abline(v=mean(thetahat.pois$alpha1), lty=2, lwd=2)
abline(h=mean(thetahat.pois$beta1), lty=2, lwd=2)
abline(v=alpha1, col="red", lty=2, lwd=2)
abline(h=beta1, col="red", lty=2, lwd=2)

## Cox
thetahat.cox <- read.csv(paste0("./sim-results/attempt", attempt, "/mcox_est.csv"))
names(thetahat.cox) <- c("iter","alpha1","beta1")

plot(beta1 ~ alpha1, data=thetahat.cox, xlim=xlimits, ylim=ylimits)
abline(v=mean(thetahat.cox$alpha1), lty=2, lwd=2)
abline(h=mean(thetahat.cox$beta1), lty=2, lwd=2)
abline(v=alpha1, col="red", lty=2, lwd=2)
abline(h=beta1, col="red", lty=2, lwd=2)


## Checking variance estimates
var.pois <- read.csv(paste0("./sim-results/attempt", attempt, "/mpois_var.csv"))
var.pois <- var.pois[, c(1, ncol(var.pois) - 1:0)]
names(var.pois) <- c("iter","alpha1","beta1")
apply(thetahat.pois[,2:3], 2, sd)
apply(sqrt(var.pois[,2:3]), 2, mean)

var.cox <- read.csv(paste0("./sim-results/attempt", attempt, "/mcox_var.csv"))
names(var.cox) <- c("iter","alpha1","beta1")
apply(thetahat.cox[,2:3], 2, sd)
apply(sqrt(var.cox[,2:3]), 2, mean)

