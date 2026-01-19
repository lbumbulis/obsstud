
source("./covariates2/source.R")
dat <- readRDS("./covariates2/test_dat_iter1.rds")

## ETM
library(etm)

# Standard approach
which(names(dat) %in% c("i","state.prev","state","u"))
names(dat)[which(names(dat) %in% c("i","state.prev","state","u"))]

etm.dat <- dat[which(dat$state.prev != dat$state), c(1,8,9,5)]
names(etm.dat) <- c("id","from","to","time")

P.logical <- (P.fn(0.5, rep(0,nvar), 0) > 0)
diag(P.logical) <- FALSE

etm.obj <- etm(
  etm.dat,
  state.names = as.character(1:12),
  tra = P.logical,
  cens.name = NULL,
  s = 0
)

tra.df <- data.frame(cbind(etm.obj$time, t(etm.obj$est[1,,]))) # take 1st slice since we start in state 1
names(tra.df) <- c("tt", paste0("state", 1:12))
tra.df <- tra.df[which(tra.df$tt <= 1),]

# Plots
grid.colour <- "grey90"

plot.AJ <- function(type) {
  age <- tra.df$tt * time.scale
  xlab <- "AGE"
  ylab <- "CUMULATIVE INCIDENCE"
  lwd <- 1.5
  
  old.par <- par(mar=c(5,4,1,1))
  
  if (type == "lung") {
    plot(age, tra.df$state7, type="s", xlab=xlab, ylab=ylab, lwd=lwd)
    lines(age, tra.df$state3, type="s", col="red", lwd=lwd)
    lines(age, tra.df$state11, type="s", col="blue", lwd=lwd)
    
    legend("topleft", col=c("red","black","blue"), lty=1, bty="n", legend=c(
      expression("(" * 0^o * ", 1)"), "(1, 1)", "(0, 1)"
    ))
  } else if (type == "lung2") {
    lung.prob <- tra.df$state3 + tra.df$state7 + tra.df$state11
    lung.prob.eversmoke <- tra.df$state7 + tra.df$state11
    
    ymax <- 0.7
    y.jump <- 0.1
    plot(age, lung.prob, type="n", xlab=xlab, ylab=ylab, ylim=c(0, ymax))
    abline(h=seq(0, ymax, by=y.jump), col=grid.colour)
    abline(v=seq(0, 120, by=10), col=grid.colour)
    lines(age, lung.prob, type="s", lwd=lwd)
    lines(age, lung.prob.eversmoke, type="s", lty=2, lwd=lwd)
    lines(age, tra.df$state11, type="s", lty=3, lwd=lwd)
    
    legend("topleft", bty="n", lty=1:3, legend=c(
      "P(Z(t) = 1), lung cancer incidence",
      expression("P(E(t)" != ~ 0^o * ", Z(t) = 1), lung cancer incidence in current or former smokers"),
      "P(E(t) = 0, Z(t) = 1), lung cancer incidence in former smokers"
    ))
  } else if (type == "lung.sum") {
    ylab <- paste0(ylab, " OF LUNG CANCER")
    lung.prob <- tra.df$state3 + tra.df$state7 + tra.df$state11
    
    plot(age, lung.prob, type="s", xlab=xlab, ylab=ylab, lwd=lwd)
  } else if (type == "lung.book") {
    plot(age, tra.book.df$state7B, type="s", xlab=xlab, ylab=ylab, lty=2, lwd=lwd)
    lines(age, tra.book.df$state7A, type="s", lwd=lwd)
    lines(age, tra.book.df$state3A, type="s", col="red", lwd=lwd)
    lines(age, tra.book.df$state3B, type="s", col="red", lty=2, lwd=lwd)
    lines(age, tra.book.df$state11A, type="s", col="blue", lwd=lwd)
    lines(age, tra.book.df$state11B, type="s", col="blue", lty=2, lwd=lwd)
    lines(age, tra.book.df$state11C, type="s", col="blue", lty=3, lwd=lwd)
    
    legend(
      "topleft", bty="n",
      col = c("red","red","black","black",rep("blue",3)),
      lty = c(1,2, 1,2, 1:3),
      legend = c(
        expression("(" * 0^o * ", 0)" ~ symbol('\256') ~  "(" * 0^o * ", 1)"),
        expression("(" * 0^o * ", 0)" ~ symbol('\256') ~  "(" * 0^o * ", 1')" ~ symbol('\256') ~  "(" * 0^o * ", 1)"),
        expression("(1, 0)" ~ symbol('\256') ~  "(1, 1)"),
        expression("(1, 0)" ~ symbol('\256') ~  "(1, 1')" ~ symbol('\256') ~  "(1, 1)"),
        expression("(0, 0)" ~ symbol('\256') ~  "(0, 1)"),
        expression("(0, 0)" ~ symbol('\256') ~  "(0, 1')" ~ symbol('\256') ~  "(0, 1)"),
        expression("(1, 0)" ~ symbol('\256') ~  "(1, 1')" ~ symbol('\256') ~  "(0, 1')" ~ symbol('\256') ~  "(0, 1)")
      )
    )
  }
  
  par(old.par)
}

plot.AJ("lung")      # CIF-lung-types.png, width=550, height=400
plot.AJ("lung2")     # CIF-lung-v2.png,    width=700, height=450 * Figure 6(a)
plot.AJ("lung.sum")  # CIF-lung-total.png, width=600, height=450
plot.AJ("lung.book") # CIF-lung-book.png,  width=600, height=450


end.states <- aggregate(state ~ i, data=dat, FUN=function(x) { x[length(x)] })
summary(factor(end.states$state))

i.absorb <- unique(end.states$i[which(!(end.states$state %in% c(1,2, 5,6, 9,10)))])

system.time(new.dat <- as.data.frame(data.table::rbindlist(lapply(i.absorb, function(i) {
  if (i %% 100 == 0) { print(i) }

  dat.i <- dat[which(dat$i==i),]
  max.ij <- max(dat.i$j)
  max.j <- J

  if (max.ij >= max.j) {
    return(NULL)
  } else {
    max.dat.i <- dat.i[max.ij,]
    extra.j <- (max.ij+1):max.j
    extra.r <- seq(from=max(dat.i$r)+1, by=1, length.out=length(extra.j))

    new.dat.i <- as.data.frame(data.table::rbindlist(rep(list(max.dat.i), length(extra.j))))
    new.dat.i$j <- extra.j
    new.dat.i$u <- extra.j/J
    new.dat.i$u.prev <- new.dat.i$u - 1/J
    new.dat.i$r <- extra.r
    return(new.dat.i)
  }
})))) # takes < 10min

dat.full <- rbind(dat, new.dat)

# Mean cumulative time spent smoking over time
mean.c <- aggregate(c ~ j, data=dat.full, FUN=mean)
idx.c <- which(mean.c$j<=J)
c.x <- mean.c$j[idx.c] * time.scale/J
c.y <- mean.c$c[idx.c] * time.scale

# Plot
old.par <- par(mar=c(5,6,1,1))

plot(
  c.x, c.y, type="n", ylim=c(0,20), xlab="AGE", # cex.axis=0.8,
  ylab=expression(atop("MEAN CUMULATIVE TIME AS", "CURRENT SMOKER (YEARS)"))
)
abline(h=seq(0, 50, by=2.5), col=grid.colour)
abline(v=seq(0, 120, by=10), col=grid.colour)
lines(c.x, c.y, type="s")

par(old.par)

# CIF-currentsmoke.png, width=500, height=450 * Figure 6(b)



