
tau <- 1
pv <- 1/2

###### SMOKING PROCESS ################

pi.w <- 0.25 # proportion of population that ever smokes

### (0^\circ) -> (1, 0)
# Baseline intensity
rho.e0.multiroot.fn <- function(logkr) {
  k <- exp(logkr[1])
  r <- exp(logkr[2])
  t.med <- (log(2))^(1/k) / r
  t.90 <- (log(10))^(1/k) / r
  # stdev <- sqrt(r^(-2) * ( gamma(1+2/k) - (gamma(1+1/k))^2 ))
  return(c(t.med - 0.2, t.90 - 0.3))
}
temp <- rootSolve::multiroot(rho.e0.multiroot.fn, start=c(0,0))
kap <- exp(temp$root[1])
rho <- exp(temp$root[2])
rho.e0.fn <- Vectorize(function(t) {
  kap*rho * (rho*t)^(kap-1)
  # same as: dweibull(t, kap, 1/rho) / pweibull(t, kap, 1/rho, lower.tail=F)
})
psi <- 5 # effect of V

### (1, 0) -> (0, 0)
lam.e00 <- 1/0.25 # baseline intensity
eta0 <- log(0.5) # effect of V

### (0, 0) -> (1, 0)
lam.e10 <- 4*lam.e00 # baseline intensity
eta1 <- log(2) # effect of V

###### FAILURE PROCESS ################
### Lung cancer
alpha1 <- log(1.25) # effect of V
beta1 <- 2*log(4)   # effect of C(t)
# Baseline intensity
lam10.multiroot.fn <- function(logkr) {
  k <- exp(logkr[1])
  r <- exp(logkr[2])
  t.med <- (log(2))^(1/k) / r
  t.90 <- (log(10))^(1/k) / r
  return(c(t.med - 0.6, t.90 - 0.8))
}
temp <- rootSolve::multiroot(lam10.multiroot.fn, start=c(0,0))
kap1 <- exp(temp$root[1])
rho1 <- exp(temp$root[2])
lam10.fn <- Vectorize(function(t) {
  kap1*rho1 * (rho1*t)^(kap1-1)
  # same as: dweibull(t, kap1, 1/rho1) / pweibull(t, kap1, 1/rho1, lower.tail=F)
})

### Death from other causes
alpha2 <- log(1.1)  # effect of V
beta2 <- 2*log(1.5) # effect of C(t)
# Baseline intensity
lam20.multiroot.fn <- function(logkr) {
  k <- exp(logkr[1])
  r <- exp(logkr[2])
  p.surv90 <- exp(-(r*0.9)^k)
  
  Lam10.fn <- function(t) { (rho1*t)^kap1 } # same as -pweibull(t, kap1, 1/rho1, lower=F, log=T)
  Lam20.fn <- function(t) { (r*t)^k }
  
  integrand0 <- Vectorize(function(s) {
    lam10.fn(s) * exp(-Lam10.fn(s) - Lam20.fn(s))
  })
  integrand1 <- Vectorize(function(s) {
    lam10.fn(s)*exp(alpha1) * exp(-Lam10.fn(s)*exp(alpha1) - Lam20.fn(s)*exp(alpha2))
  })
  pi.c <- (1-pv)*integrate(integrand0, 0, tau)$value + pv*integrate(integrand1, 0, tau)$value
  
  return(c(p.surv90 - 0.1, pi.c - 0.1))
}
temp <- rootSolve::multiroot(lam20.multiroot.fn, start=c(0,0))
kap2 <- exp(temp$root[1])
rho2 <- exp(temp$root[2])
# lam20.fn <- Vectorize(function(t) {
#   kap2*rho2 * (rho2*t)^(kap2-1)
#   # same as: dweibull(t, kap2, 1/rho2) / pweibull(t, kap2, 1/rho2, lower.tail=F)
# })
lam20.fn <- Vectorize(function(t) { 0 }) # TODO: switch to nonzero

###### TRANSITION MATRICES ############

# One-step transition probability matrix, given time tt, cumulative exposure cc, and covariates v and w
P.fn <- function(tt, cc, v, w) {
  P <- matrix(c(
    0, w*rho.e0.fn(tt)*exp(psi*v), 0, lam10.fn(tt)*exp(alpha1*v), 0, 0,
    0, 0, lam.e00*exp(eta0*v), 0, lam10.fn(tt)*exp(alpha1*v + beta1*cc), 0,
    0, lam.e10*exp(eta1*v), 0, 0, 0, lam10.fn(tt)*exp(alpha1*v + beta1*cc),
    rep(0, 6*3)
  ), nrow=2*3, byrow=T) * 1/R
  diag(P) <- 1 - rowSums(P)
  
  return(P)
}

###### SIMULATION PARAMETERS ##########

R <- 1000  # size of partition of [0,tau)
n <- 1000 # sample size

nsim <- 1000 # number of simulation replicates

###### DATA GENERATION ################

get.E.for.state <- function(state) {
  # -1 represents 0^\circ
  ifelse(state %in% c(1,4), -1, ifelse(state %in% c(2,5), 1, 0))
}
get.N.for.state <- function(state) {
  as.numeric(state >= 4)
}

generate.data <- function() {
  u <- (1:R)/R
  v.full <- rbinom(n, size=1, prob=pv)
  w.full <- rbinom(n, size=1, prob=pi.w)
  
  dat <- data.frame(
    i = rep(1:n, each=R),
    v = rep(v.full, each=R),
    w = rep(w.full, each=R),
    r = rep(1:R, times=n),
    u = rep(u, times=n),
    E = NA,
    N.bar = NA,
    state.prev = NA,
    state = NA,
    c = NA
  )
  
  P.onset00 <- t(sapply(1:R, function(r) {
    P.fn(u[r], cc=0, v=0, w=0)[1,]
  }))
  P.onset01 <- t(sapply(1:R, function(r) {
    P.fn(u[r], cc=0, v=0, w=1)[1,]
  }))
  P.onset10 <- t(sapply(1:R, function(r) {
    P.fn(u[r], cc=0, v=1, w=0)[1,]
  }))
  P.onset11 <- t(sapply(1:R, function(r) {
    P.fn(u[r], cc=0, v=1, w=1)[1,]
  }))
  
  system.time(sapply(1:n, function(i) {
    print(paste0(Sys.time(), ": ", i))
    idx <- which(dat$i==i)
    v <- dat$v[idx[1]]
    w <- dat$w[idx[1]]
    
    # Multinomial draws via inverse CDF method, to determine 1st time when not in (0^\circ, 0)
    unif <- runif(R)
    if (v==0 && w==0) {
      P <- P.onset00
    } else if (v==0 && w==1) {
      P <- P.onset01
    } else if (v==1 && w==0) {
      P <- P.onset10
    } else if (v==1 && w==1) {
      P <- P.onset11
    }
    cumP <- t(apply(P, 1, cumsum))
    onset.states <- max.col(cumP >= unif, ties.method="first")
    
    if (any(onset.states!=1)) {
      onset.idx <- min(which(onset.states!=1)) # 1st time when not in (0^\circ, 0)
      onset.state <- onset.states[onset.idx]   # state occupied at that time
      
      pre.onset.idx <- idx[1:(onset.idx-1)]
      dat$E[pre.onset.idx] <<- -1
      dat$N.bar[pre.onset.idx] <<- 0
      dat$state.prev[pre.onset.idx] <<- 1
      dat$state[pre.onset.idx] <<- 1
      dat$c[pre.onset.idx] <<- 0
      
      E.prev <- get.E.for.state(onset.state)
      N.bar <- get.N.for.state(onset.state)
      c.prev <- as.numeric(E.prev == 1)/R
      state.prev <- onset.state # 1=(0^\circ,0), 2=(1,0), 3=(0,0), 4=(0^\circ,1), 5=(1,1), 6=(0,1)
      
      at.onset.idx <- idx[onset.idx]
      dat$E[at.onset.idx] <<- E.prev
      dat$N.bar[at.onset.idx] <<- N.bar
      dat$state.prev[at.onset.idx] <<- 1
      dat$state[at.onset.idx] <<- state.prev
      dat$c[at.onset.idx] <<- c.prev
      
      r <- onset.idx + 1
      while (N.bar==0 && r <= R) {
        P <- P.fn(u[r], c.prev, v, w)
        state <- sample(ncol(P), size=1, prob=P[state.prev,])
        E <- ifelse(state %in% c(1,4), -1, ifelse(state %in% c(2,5), 1, 0))
        N.bar <- as.numeric(state >= 4)
        c <- c.prev + as.numeric(E == 1)/R
        
        dat$E[idx[r]] <<- E
        dat$N.bar[idx[r]] <<- N.bar
        dat$state.prev[idx[r]] <<- state.prev
        dat$state[idx[r]] <<- state
        dat$c[idx[r]] <<- c
        
        E.prev <- E
        c.prev <- c
        state.prev <- state
        r <- r+1
      }
      
      if (exists("E")) {
        post.fail.idx <- idx[r:R]
      } else if (onset.idx < R) {
        post.fail.idx <- idx[(onset.idx+1):R]
      }
      if (exists("post.fail.idx")) {
        dat$E[post.fail.idx] <<- ifelse(exists("E"), E, -1)
        dat$N.bar[post.fail.idx] <<- 1
        dat$state.prev[post.fail.idx] <<- ifelse(exists("E"), state, 4)
        dat$state[post.fail.idx] <<- ifelse(exists("E"), state, 4)
        dat$c[post.fail.idx] <<- ifelse(exists("E"), c, 0)
      }
      
    } else {
      dat$E[idx] <<- -1
      dat$N.bar[idx] <<- 0
      dat$state.prev[idx] <<- 1
      dat$state[idx] <<- 1
      dat$c[idx] <<- 0
    }
  }))
  
  # dat <- dat[which(!is.na(dat$E)),]
  dat$u.prev <- dat$u - 1/R
  
  return(dat)
}


