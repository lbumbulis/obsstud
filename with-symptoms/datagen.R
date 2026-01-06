
# One-step transition intensity matrix, given time tt,
# cumulative exposure variables x, and covariate v
P.fn <- function(tt, x, v) {
  P <- matrix(c(
    ## e = 0^\circ
    # z = 0
    0,
    lam1p.fn(tt)*exp(t(x) %*% beta1p + v*gamma1p),
    lam1.fn(tt)*exp(t(x) %*% beta1 + v*gamma1),
    lam2.fn(tt)*exp(t(x) %*% beta2 + v*gamma2),
    psi1c.fn(tt)*exp(v*eta1c), 0, 0, 0,
    rep(0, 4),
    # z = 1'
    0, 0, lam1.fn(tt)*exp(t(x) %*% beta1 + v*gamma1), 0,
    rep(0, 4*2),
    # z = 1 or 2
    rep(0, 4*3*2),
    
    ## e = 1
    # z = 0
    rep(0, 4),
    0,
    lam1p.fn(tt)*exp(t(x) %*% beta1p + v*gamma1p),
    lam1.fn(tt)*exp(t(x) %*% beta1 + v*gamma1),
    lam2.fn(tt)*exp(t(x) %*% beta2 + v*gamma2),
    psi0*exp(t(x) %*% alpha0 + v*eta0), 0, 0, 0,
    # z = 1'
    rep(0, 4),
    0, 0, lam1.fn(tt)*exp(t(x) %*% beta1 + v*gamma1), 0,
    0, psi0*exp(t(x) %*% alpha0 + v*eta0), 0, 0,
    # z = 1 or 2
    rep(0, 4*3*2),
    
    ## e = 0
    # z = 0
    rep(0, 4),
    psi1*exp(t(x) %*% alpha1 + v*eta1), 0, 0, 0,
    0,
    lam1p.fn(tt)*exp(t(x) %*% beta1p + v*gamma1p),
    lam1.fn(tt)*exp(t(x) %*% beta1 + v*gamma1),
    lam2.fn(tt)*exp(t(x) %*% beta2 + v*gamma2),
    # z = 1'
    rep(0, 4*2),
    0, 0, lam1.fn(tt)*exp(t(x) %*% beta1 + v*gamma1), 0,
    # z = 1 or 2
    rep(0, 4*3*2)
    
  ), nrow=4*3, byrow=T) * 1/R
  diag(P) <- 1 - rowSums(P)
  
  return(P)
}

get.E.for.state <- function(state) {
  # -1 represents 0^\circ
  ifelse(state %in% 1:4, -1, ifelse(state %in% 5:8, 1, 0))
}
get.Z.for.state <- function(state) {
  # 10 represents 1'
  ifelse(state %% 4 == 2, 10, ifelse(state %% 4 == 1, 0, (state-2) %% 4))
}

get.exit.info <- function(tt, cc, b.start, v, state) {
  surv.prob <- runif(1)
  
  integrand.1p <- Vectorize(function(u) { lam1p.fn(u) * exp(beta1p[2]*u )})
  integrand.1 <- Vectorize(function(u) { lam1.fn(u) * exp(beta1[2]*u) })
  integrand.2 <- Vectorize(function(u) { lam2.fn(u) * exp(beta2[2]*u) })
  
  get.x4.idx <- function(u) {
    b <- b.start + u - tt # time since most recent quit
    return(findInterval(b, b.breaks))
  }
  
  integrand.1p.b <- Vectorize(function(u) { lam1p.fn(u) * exp(beta1p[3+get.x4.idx(u)]) })
  integrand.1.b <- Vectorize(function(u) { lam1.fn(u) * exp(beta1[3+get.x4.idx(u)]) })
  integrand.2.b <- Vectorize(function(u) { lam2.fn(u) * exp(beta2[3+get.x4.idx(u)]) })
  
  E <- get.E.for.state(state)
  Z <- get.Z.for.state(state)
  
  # Determine exit time
  if (E==-1 && Z==0) {
    uniroot.fn <- function(s) {
      exp(-(
        exp(gamma1p*v) * integrate(lam1p.fn, lower=tt, upper=tt+s)$value +
          exp(gamma1*v) * integrate(lam1.fn, lower=tt, upper=tt+s)$value +
          exp(gamma2*v) * integrate(lam2.fn, lower=tt, upper=tt+s)$value +
          exp(eta1c*v) * integrate(psi1c.fn, lower=tt, upper=tt+s)$value
      )) - surv.prob
    }
  } else if (E==-1 && Z==10) {
    uniroot.fn <- function(s) {
      exp(-(
        exp(beta1[nvar] + gamma1*v) * integrate(lam1.fn, lower=tt, upper=tt+s)$value
      )) - surv.prob
    }
  } else if (E==1 && Z==0) {
    uniroot.fn <- function(s) {
      exp(-(
        exp(beta1p[1] + (cc-tt)*beta1p[2] + gamma1p*v) *
            integrate(integrand.1p, lower=tt, upper=tt+s)$value +
          exp(beta1[1] + (cc-tt)*beta1[2] + gamma1*v) *
            integrate(integrand.1, lower=tt, upper=tt+s)$value +
          exp(beta2[1] + (cc-tt)*beta2[2] + gamma2*v) *
            integrate(integrand.2, lower=tt, upper=tt+s)$value +
          s * psi0*exp(eta0*v)
      )) - surv.prob
    }
  } else if (E==1 && Z==10) {
    uniroot.fn <- function(s) {
      exp(-(
        exp(beta1[1] + (cc-tt)*beta1[2] + beta1[nvar] + gamma1*v) *
            integrate(integrand.1, lower=tt, upper=tt+s)$value +
          s * psi0*exp(alpha0[nvar] + eta0*v)
      )) - surv.prob
    }
  } else if (E==0 && Z==0) {
    uniroot.fn <- function(s) {
      exp(-(
        exp(beta1p[3] + gamma1p*v) * integrate(integrand.1p.b, lower=tt, upper=tt+s)$value +
          exp(beta1[3] + gamma1*v) * integrate(integrand.1.b, lower=tt, upper=tt+s)$value +
          exp(beta2[3] + gamma2*v) * integrate(integrand.2.b, lower=tt, upper=tt+s)$value +
          s * psi1*exp(alpha1[3] + eta1*v)
      )) - surv.prob
    }
  } else if (E==0 && Z==10) {
    uniroot.fn <- function(s) {
      exp(-(
        exp(beta1[3] + beta1[nvar] + gamma1*v) * integrate(integrand.1.b, lower=tt, upper=tt+s)$value
      )) - surv.prob
    }
  }
  
  if (uniroot.fn(1/R) < 0) {
    exit.time <- tt + 1/R
  } else {
    exit.time <- tt + uniroot(uniroot.fn, interval=c(1/R, 2))$root
  }
  
  if (exit.time > tau) {
    is.exit <- F
    exit.idx <- exit.state <- NA
  } else {
    is.exit <- T
    
    # Determine exit state
    x5 <- as.numeric(Z==10)
    if (E==-1) {
      x <- c(rep(0, nvar-1), x5)
    } else if (E==1) {
      x <- c(1, cc + exit.time - tt, rep(0, nvar-3), x5)
    } else if (E==0) {
      x4 <- rep(0, RB)
      x4[get.x4.idx(exit.time)] <- 1
      x <- c(0, 0, 1, x4, x5)
    }
    
    P <- P.fn(exit.time, x, v)
    diag(P) <- 0 # prevent "transitions" where they stay in the same state
    exit.state <- which(rmultinom(n=1, size=1, prob=P[state,]) > 0)
    
    # Discretize exit time
    exit.idx <- ceiling(exit.time * R)
  }
  
  return(list(is.exit=is.exit, exit.idx=exit.idx, exit.state=exit.state))
}

generate.data <- function(n=nn) {
  v.full <- rbinom(n, size=1, prob=pv)
  
  dat <- data.frame(
    i = rep(1:n, each=R),
    v = rep(v.full, each=R),
    r = rep(1:R, times=n),
    u = rep((1:R)/R, times=n),
    E = NA,
    Z = NA,
    state.prev = NA,
    state = NA,
    c = NA,
    b = NA
  )
  
  fill.dat <- function(idx, E, Z, state.prev, state, cc, b) {
    dat$E[idx] <<- E
    dat$Z[idx] <<- Z
    dat$state.prev[idx] <<- state.prev
    dat$state[idx] <<- state
    dat$c[idx] <<- cc
    dat$b[idx] <<- b
  }
  
  system.time(sapply(1:n, function(i) {
    print(paste0(Sys.time(), ": ", i))
    i.idx <- which(dat$i==i)
    v <- dat$v[i.idx[1]]
    
    state <- 1
    E <- get.E.for.state(state)
    Z <- get.Z.for.state(state)
    cc <- 0
    b <- Inf
    r <- 1
    
    while (Z %in% c(0,10) && r <= R) { # Z=10 means Z=1'
      exit.info <- get.exit.info((r-1)/R, cc, b, v, state)
      exit.idx <- exit.info$exit.idx
      
      if (!exit.info$is.exit) {
        break
      } else {
        exit.state <- exit.info$exit.state
        
        if (E == 1) { # currently smoking
          if (exit.idx > r) {
            fill.dat(
              i.idx[r:(exit.idx-1)],
              E=E, Z=Z, state.prev=state, state=state,
              c = cc + seq(1/R, by=1/R, length.out=exit.idx-r),
              b = 0
            )
          }
          cc <- cc + (exit.idx - r)/R
          b <- 0
        } else if (E == 0) { # previously smoked
          if (exit.idx > r) {
            fill.dat(
              i.idx[r:(exit.idx-1)],
              E=E, Z=Z, state.prev=state, state=state, c=cc,
              b = b + seq(1/R, by=1/R, length.out=exit.idx-r)
            )
          }
          b <- b + (exit.idx - r)/R
        } else { # (E == -1) never smoked
          if (exit.idx > r) {
            fill.dat(
              i.idx[r:(exit.idx-1)],
              E=E, Z=Z, state.prev=state, state=state, c=cc, b=b # c=b=0 here
            )
          }
        }
        
        E <- get.E.for.state(exit.state)
        Z <- get.Z.for.state(exit.state)
        fill.dat(
          i.idx[exit.idx],
          E = E,
          Z = Z,
          state.prev = state,
          state = exit.state,
          c = cc,
          b = b
        )
        state <- exit.state
        r <- exit.idx + 1
      }
    }
    
    # Current state is never exited; fill in remaining times accordingly
    if (r <= R) {
      c.vec <- cc
      b.vec <- b
      
      if (Z %in% c(0,10)) { # unfailed
        if (E==1) {
          c.vec <- c.vec + seq(1/R, by=1/R, length.out=R-r+1)
        } else if (E==0) {
          b.vec <- b.vec + seq(1/R, by=1/R, length.out=R-r+1)
        }
      }
      
      fill.dat(
        i.idx[r:R],
        E = E,
        Z = Z,
        state.prev = state,
        state = state,
        c = c.vec,
        b = b.vec
      )
    }
  }))
  
  dat$u.prev <- dat$u - 1/R
  
  return(dat)
}
