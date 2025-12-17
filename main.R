
library(ggplot2)

source("source.R")

attempt <- 3

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

plot.dat <- as.data.frame(Reduce("+", all.state.occupancy) / length(all.state.occupancy) / R)
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
hist(quit.dat, xlab="mean # quit attempts", main="")
mean(quit.dat) # 0.912

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

# Mean pre-failure smoking duration among ever-smokers
dur.dat <- sapply(1:length(all.smoke.dur), function(iter) {
  temp <- all.smoke.dur[[iter]]
  mean(temp$c[which(temp$c > 0)])
})
hist(dur.dat, xlab="mean smoking duration", main="")
mean(dur.dat) # 0.330

