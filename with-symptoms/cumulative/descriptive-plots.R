
# Testing data generation
source("./with-symptoms/source.R")
source("./with-symptoms/cumulative/datagen.R")

sim.seeds <- readRDS("sim_seeds_nsim1000.rds")

#######################################
n <- nn
.Random.seed <- sim.seeds[[1]]
system.time(dat <- generate.data(n))

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

# Bookkeeping
dat$alive.smoke <- as.numeric(dat$E == 1 & dat$Z %in% c(0,10))
smoke.prob <- aggregate(alive.smoke ~ j, data=dat, FUN=mean)
sum(smoke.prob$alive.smoke) * time.scale/J
# Mean time spent smoking prior to failure/censoring is ~21.2 years

dat$alive.quit <- as.numeric(dat$E == 0 & dat$Z %in% c(0,10))
quit.prob <- aggregate(alive.quit ~ j, data=dat, FUN=mean)
sum(smoke.prob$alive.quit) * time.scale/J / 365
# Mean time spent as former smoker prior to failure/censoring is ~0 days

dat$is.symptom.quit <- as.numeric(dat$state.prev==6 & dat$state==10)
ever.symptom.quit <- aggregate(is.symptom.quit ~ i, data=dat, FUN=max)
names(ever.symptom.quit)[2] <- "ever.symptom.quit"
dat <- merge(dat, ever.symptom.quit)

dat$state.book <- dat$state
dat$state.book[which(dat$state==3 & dat$state.prev==1)] <- "3A"
dat$state.book[which(dat$state==3 & dat$state.prev==2)] <- "3B"
dat$state.book[which(dat$state==7 & dat$state.prev==5)] <- "7A"
dat$state.book[which(dat$state==7 & dat$state.prev==6)] <- "7B"
dat$state.book[which(dat$state==11 & dat$state.prev==9)] <- "11A"
dat$state.book[which(
  dat$state==11 & dat$state.prev==10 & dat$ever.symptom.quit==0
)] <- "11B"
dat$state.book[which(
  dat$state==11 & dat$state.prev==10 & dat$ever.symptom.quit==1
)] <- "11C"

which(names(dat)=="state.book")

etm.book <- dat[which(dat$state.prev != dat$state), c(1,8,27,5)]
names(etm.book) <- c("id","from","to","time")
etm.book$from <- as.character(etm.book$from)

book.states <- c("1","2","3A","3B","4", "5","6","7A","7B","8", "9","10","11A","11B","11C","12")

etm.book.obj <- etm(
  etm.book,
  state.names = book.states,
  tra = matrix(as.logical(c(
    # exiting e=0^\circ
    0,1,1,0,1, 1,0,0,0,0, 0,0,0,0,0,0,
    0,0,0,1,0, 0,0,0,0,0, 0,0,0,0,0,0,
    rep(0, 3*16),
    # exiting e=1
    0,0,0,0,0, 0,1,1,0,1, 1,0,0,0,0,0,
    0,0,0,0,0, 0,0,0,1,0, 0,1,0,0,0,0,
    rep(0, 3*16),
    # exiting e=0
    0,0,0,0,0, 1,0,0,0,0, 0,1,1,0,0,1,
    0,0,0,0,0, 0,0,0,0,0, 0,0,0,1,1,0,
    rep(0, 4*16)
  )), nrow=12+4, byrow=T),
  cens.name = NULL,
  s = 0
)

tra.book.df <- data.frame(cbind(etm.book.obj$time, t(etm.book.obj$est[1,,])))
names(tra.book.df) <- c("tt", paste0("state", book.states))

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
    
    plot(age, lung.prob, type="n", xlab=xlab, ylab=ylab, ylim=c(0, 0.2))
    abline(h=seq(0, 0.2, by=0.025), col=grid.colour)
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
plot.AJ("lung2")     # CIF-lung.png,       width=700, height=450 *
plot.AJ("lung.sum")  # CIF-lung-total.png, width=600, height=450
plot.AJ("lung.book") # CIF-lung-book.png,  width=600, height=450




## Mean cumulative time spent in quit state over time
dat$is.quit <- (dat$state.prev==5 & dat$state==9) | (dat$state.prev==6 & dat$state==10)
foo <- aggregate(is.quit ~ i, data=dat, FUN=sum)
max(foo$is.quit) # no one individual quits more than once

quitters <- unique(dat$i[which(is.finite(dat$b) & dat$b>0)])
dat$max.b <- dat$b
for (i in quitters) {
  i.idx <- which(dat$i==i)
  dat.i <- dat[i.idx,]
  max.b <- max(dat.i$b[which(is.finite(dat.i$b))])
  max.b.idx <- which(dat.i$b==max.b)
  
  if (max(max.b.idx) != length(i.idx)) { print(i) }
}
# Nobody resumes smoking

# Confirm:
sum(dat$state.prev==9 & dat$state==5) # 0

# So we can approach cumulative time in quit state exactly as with time in the
# smoking state ... but first we need to pad the end of the dataframe with time
# past follow-up for those who end in an absorbing state.
end.states <- aggregate(state ~ i, data=dat, FUN=function(x) { x[length(x)] })
summary(factor(end.states$state))

i.absorb <- unique(end.states$i[which(!(end.states$state %in% c(1,2, 5,6, 9,10)))])

system.time(new.dat <- as.data.frame(data.table::rbindlist(lapply(i.absorb, function(i) {
  if (i %% 100 == 0) { print(i) }
  
  dat.i <- dat[which(dat$i==i),]
  max.ij <- max(dat.i$j)
  max.j <- (A+tau)*J
  
  if (max.ij==max.j) {
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
})))) # takes about 5min

dat.full <- rbind(dat, new.dat)

## Now make the plot
dat.full$b.finite <- dat.full$b
dat.full$b.finite[which(!is.finite(dat.full$b))] <- 0
mean.b <- aggregate(b.finite ~ j, data=dat.full, FUN=mean)
idx.b <- which(mean.b$j<=850)
b.x <- mean.b$j[idx.b] * time.scale/J
b.y <- mean.b$b.finite[idx.b] * time.scale * 365

# Mean cumulative time spent smoking over time
mean.c <- aggregate(c ~ j, data=dat.full, FUN=mean)
idx.c <- which(mean.c$j<=1000)
c.x <- mean.c$j[idx.c] * time.scale/J
c.y <- mean.c$c[idx.c] * time.scale

# Plot
old.par <- par(mfrow=c(1,2), mar=c(5,6,1,1))

plot(
  c.x, c.y, type="n", cex.axis=0.8, xlab="AGE",
  ylab=expression(atop("MEAN CUMULATIVE TIME AS", "CURRENT SMOKER (YEARS)"))
)
abline(h=seq(0, 20, by=2.5), col=grid.colour)
abline(v=seq(0, 120, by=10), col=grid.colour)
lines(c.x, c.y, type="s")

plot(
  b.x, b.y, type="n", cex.axis=0.8, xlab="AGE",
  ylab=expression(atop("MEAN CUMULATIVE TIME AS", "FORMER SMOKER (DAYS)"))
)
abline(h=seq(0, 120, by=10), col=grid.colour)
abline(v=seq(0, 120, by=10), col=grid.colour)
lines(b.x, b.y, type="s")

par(old.par)

# CIF-smoke.png, width=750, height=350



