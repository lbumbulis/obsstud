
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

# P is a matrix with one column for each state and one row for each time point
# (i.e., vector of multinomial probabilities)
my.rmultinom <- function(P) {
  unif <- runif(nrow(P))
  cumP <- t(apply(P, 1, cumsum))
  return(max.col(cumP >= unif, ties.method="first"))
}

get.exit.info <- function(start.r, c.prev, v, state.prev) {
  if (state.prev==2) {
    P.exit <- t(sapply(start.r:R, function(rr) {
      P.fn(rr/R, cc=c.prev+(rr-start.r+1)/R, v=v, w=1)[state.prev,]
    }))
  } else if (state.prev==3) {
    P.exit <- t(sapply(start.r:R, function(rr) {
      P.fn(rr/R, cc=c.prev, v=v, w=1)[state.prev,]
    }))
  }
  exit.states <- my.rmultinom(P.exit)
  
  if (any(exit.states != state.prev)) {
    is.exit <- T
    exit.idx <- min(which(exit.states != state.prev)) # next time when not in current state
    exit.state <- exit.states[exit.idx]               # state occupied at that time
  } else {
    is.exit <- F
    exit.idx <- exit.state <- NA
  }
  return(list(is.exit=is.exit, exit.idx=exit.idx, exit.state=exit.state))
}

get.exit.info.cts <- function(start.r, c.prev, v, state.prev) {
  surv.prob <- runif(1)
  if (state.prev == 2) {
    lam.nonfail <- lam.e00 * exp(eta0*v) # quit intensity
  } else if (state.prev == 3) {
    lam.nonfail <- lam.e10 * exp(eta1*v) # relapse intensity
  }
  
  # Determine exit time
  if (state.prev == 2) {
    integrate.fn <- Vectorize(function(u) {
      lam10.fn(u) * exp(beta1*u)
    })
    uniroot.fn <- function(s) {
      exp(
        -integrate(integrate.fn, lower=start.r/R, upper=start.r/R+s)$value *
          exp(alpha1*v + beta1*(c.prev-start.r/R)) -
          lam.nonfail * s
      ) - surv.prob
    }
  } else if (state.prev == 3) {
    uniroot.fn <- function(s) {
      exp(
        -(integrate(lam10.fn, lower=start.r/R, upper=start.r/R+s)$value *
            exp(alpha1*v + beta1*c.prev)) -
          lam.nonfail * s
      ) - surv.prob
    }
  }
  exit.time <- start.r/R + uniroot(uniroot.fn, interval=c(0, 2))$root
  
  if (exit.time > tau) {
    is.exit <- F
    exit.idx <- exit.state <- NA
  } else {
    is.exit <- T
    
    # Determine exit state
    c.exit <- c.prev
    if (state.prev == 2) {
      c.exit <- c.exit + exit.time - start.r/R
    }
    is.fail <- (rbinom(
      1, size=1,
      prob = (
        lam.nonfail /
          (lam.nonfail + lam10.fn(exit.time)*exp(alpha1*v + beta1*c.exit))
      )
    ) == 0)
    if (is.fail) {
      exit.state <- ifelse(state.prev==2, 5, 6)
    } else {
      exit.state <- ifelse(state.prev==2, 3, 2)
    }
    
    # Discretize exit time
    exit.idx <- ceiling(exit.time * R) - start.r
  }
  
  return(list(is.exit=is.exit, exit.idx=exit.idx, exit.state=exit.state))
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
  
  P.onset00 <- t(sapply(1:R, function(rr) { P.fn(u[rr], cc=0, v=0, w=0)[1,] }))
  P.onset01 <- t(sapply(1:R, function(rr) { P.fn(u[rr], cc=0, v=0, w=1)[1,] }))
  P.onset10 <- t(sapply(1:R, function(rr) { P.fn(u[rr], cc=0, v=1, w=0)[1,] }))
  P.onset11 <- t(sapply(1:R, function(rr) { P.fn(u[rr], cc=0, v=1, w=1)[1,] }))
  
  fill.dat <- function(idxx, E, N.bar, state.prev, state, cc) {
    dat$E[idxx] <<- E
    dat$N.bar[idxx] <<- N.bar
    dat$state.prev[idxx] <<- state.prev
    dat$state[idxx] <<- state
    dat$c[idxx] <<- cc
  }
  
  system.time(sapply(1:n, function(i) {
    # print(paste0(Sys.time(), ": ", i))
    idx <- which(dat$i==i)
    v <- dat$v[idx[1]]
    w <- dat$w[idx[1]]
    
    # Multinomial draws via inverse CDF method, to determine 1st time when not in (0^\circ, 0)
    if (v==0 && w==0) {
      P <- P.onset00
    } else if (v==0 && w==1) {
      P <- P.onset01
    } else if (v==1 && w==0) {
      P <- P.onset10
    } else if (v==1 && w==1) {
      P <- P.onset11
    }
    onset.states <- my.rmultinom(P) # 1=(0^\circ,0), 2=(1,0), 3=(0,0), 4=(0^\circ,1), 5=(1,1), 6=(0,1)
    
    # If the person ever leaves state (0^\circ, 0) ...
    if (any(onset.states!=1)) {
      onset.idx <- min(which(onset.states!=1)) # 1st time when not in (0^\circ, 0)
      onset.state <- onset.states[onset.idx]   # state occupied at that time
      
      if (onset.idx > 1) {
        fill.dat(idx[1:(onset.idx-1)], E=-1, N.bar=0, state.prev=1, state=1, cc=0)
      }
      
      E.prev <- get.E.for.state(onset.state)
      N.bar <- get.N.for.state(onset.state)
      c.prev <- 0
      state.prev <- onset.state
      
      fill.dat(idx[onset.idx], E=E.prev, N.bar=N.bar, state.prev=1, state=state.prev, cc=c.prev)
      
      r <- onset.idx + 1
      while (N.bar==0 && r <= R) {
        exit.info <- get.exit.info.cts(r, c.prev, v, state.prev)
        exit.idx <- exit.info$exit.idx
        
        if (!exit.info$is.exit) {
          break
        } else {
          exit.state <- exit.info$exit.state
          
          # If previous state is (1,0) ...
          if (state.prev==2) {
            if (exit.idx > 1) {
              fill.dat(
                idx[r:(r+exit.idx-2)],
                E=1, N.bar=0, state.prev=2, state=2,
                c = c.prev + seq(1/R, by=1/R, length.out=exit.idx-1)
              )
            }
            
            E.prev <- get.E.for.state(exit.state)
            N.bar <- get.N.for.state(exit.state)
            c.prev <- c.prev + exit.idx/R
            state.prev <- exit.state
            
            fill.dat(
              idx[r+exit.idx-1],
              E = E.prev,
              N.bar = N.bar,
              state.prev = 2,
              state = state.prev,
              c = c.prev
            )
            
            r <- r + exit.idx
            # ... else if previous state is (0,0) ...
          } else if (state.prev==3) {
            if (exit.idx > 1) {
              fill.dat(idx[r:(r+exit.idx-2)], E=0, N.bar=0, state.prev=3, state=3, c=c.prev)
            }
            
            E.prev <- get.E.for.state(exit.state)
            N.bar <- get.N.for.state(exit.state)
            state.prev <- exit.state
            
            fill.dat(
              idx[r+exit.idx-1],
              E = E.prev,
              N.bar = N.bar,
              state.prev = 3,
              state = state.prev,
              c = c.prev
            )
            
            r <- r + exit.idx
          } else { stop() }
        }
      }
      
      # Current state is never exited; fill in remaining times accordingly
      if (r <= R) {
        if (E.prev==1 && N.bar==0) {
          c.vec <- c.prev + seq(1/R, by=1/R, length.out=R-r+1)
        } else {
          c.vec <- c.prev
        }
        fill.dat(
          idx[r:R],
          E = E.prev,
          N.bar = N.bar,
          state.prev = state.prev,
          state = state.prev,
          cc = c.vec
        )
      }
      
    } else {
      fill.dat(idx, E=-1, N.bar=0, state.prev=1, state=1, cc=0)
    }
  }))
  
  dat$u.prev <- dat$u - 1/R
  
  return(dat)
}


