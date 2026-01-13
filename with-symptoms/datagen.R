
# One-step transition intensity matrix for entering time tt (1 unit = 100 years),
# given cumulative exposure variables x and covariate v
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
    
  ), nrow=4*3, byrow=T) * 1/J # times 1/J since that is the discretization interval length
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

# P is a matrix with one column for each state and one row for each time point
# (i.e., vector of multinomial probabilities)
my.rmultinom <- function(P) {
  unif <- runif(nrow(P))
  cumP <- t(apply(P, 1, cumsum))
  return(max.col(cumP >= unif, ties.method="first"))
}

get.sojourn.time <- function(uniroot.fn) {
  if (uniroot.fn(1e-10) < 0) {
    sojourn.time <- 1e-10
  } else if (uniroot.fn(1) > 0) {
    sojourn.time <- 1
  } else {
    sojourn.time <- uniroot(uniroot.fn, interval=c(1e-10, 1))$root
  }
  return(sojourn.time)
}

get.first.exit.time <- Vectorize(function(v, surv.prob) {
  uniroot.fn <- function(s) {
    exp(-(
      exp(gamma1p*v) * integrate(lam1p.fn, lower=0, upper=s)$value +
        exp(gamma1*v) * integrate(lam1.fn, lower=0, upper=s)$value +
        exp(gamma2*v) * integrate(lam2.fn, lower=0, upper=s)$value +
        exp(eta1c*v) * integrate(psi1c.fn, lower=0, upper=s)$value
    )) - surv.prob
  }
  
  return(get.sojourn.time(uniroot.fn))
})

get.first.exit.state <- Vectorize(function(exit.time, v) {
  P <- P.fn(exit.time, x=rep(0, nvar), v)
  diag(P) <- 0 # prevent "transitions" where they stay in the same state
  return(which(rmultinom(n=1, size=1, prob=P[1,]) > 0))
})

# This function applies to transitions after the first transition out of (0^\circ, 0)
get.exit.info <- function(j, max.idx, cc, b.start, v, state, method) {
  E <- get.E.for.state(state)
  Z <- get.Z.for.state(state)
  
  x6 <- as.numeric(Z==10)
  if (E==-1) {
    x <- c(rep(0, nvar-1), x6)
  } else if (E==1) {
    x <- c(1, NA, rep(0, nvar-3), x6)
  } else if (E==0 && Z==0) {
    x <- c(0, 0, 1, rep(NA, RB), 0, x6)
  } else if (E==0 && Z==10) {
    x <- c(0, 0, 1, rep(0, RB), cc*time.scale, x6)
  }
  
  if (method == "cts") {
    surv.prob <- runif(1)
    
    integrand.1p <- Vectorize(function(u) { lam1p.fn(u) * exp(beta1p[2]*u) })
    integrand.1 <- Vectorize(function(u) { lam1.fn(u) * exp(beta1[2]*u) })
    integrand.2 <- Vectorize(function(u) { lam2.fn(u) * exp(beta2[2]*u) })
    
    tt <- (j-1)/J
    get.x4.idx <- function(u) {
      b <- b.start + u - tt # time since most recent quit
      return(findInterval(b, b.breaks))
    }
    
    integrand.1p.b <- Vectorize(function(u) { lam1p.fn(u) * exp(beta1p[3+get.x4.idx(u)]*cc*time.scale) })
    integrand.1.b <- Vectorize(function(u) { lam1.fn(u) * exp(beta1[3+get.x4.idx(u)]*cc*time.scale) })
    integrand.2.b <- Vectorize(function(u) { lam2.fn(u) * exp(beta2[3+get.x4.idx(u)]*cc*time.scale) })
    
    # Determine exit time
    if (E==-1 && Z==10) {
      uniroot.fn <- function(s) {
        exp(-(
          exp(beta1[nvar] + gamma1*v) * integrate(lam1.fn, lower=tt, upper=tt+s)$value
        )) - surv.prob
      }
    } else if (E==1 && Z==0) {
      uniroot.fn <- function(s) {
        exp(-(
          exp(beta1p[1] + (cc-tt)*time.scale*beta1p[2] + gamma1p*v) *
            integrate(integrand.1p, lower=tt, upper=tt+s)$value +
            exp(beta1[1] + (cc-tt)*time.scale*beta1[2] + gamma1*v) *
            integrate(integrand.1, lower=tt, upper=tt+s)$value +
            exp(beta2[1] + (cc-tt)*time.scale*beta2[2] + gamma2*v) *
            integrate(integrand.2, lower=tt, upper=tt+s)$value +
            s * psi0*exp(eta0*v)
        )) - surv.prob
      }
    } else if (E==1 && Z==10) {
      uniroot.fn <- function(s) {
        exp(-(
          exp(beta1[1] + (cc-tt)*time.scale*beta1[2] + beta1[nvar] + gamma1*v) *
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
            s * psi1*exp(eta1*v)
        )) - surv.prob
      }
    } else if (E==0 && Z==10) {
      uniroot.fn <- function(s) {
        exp(-(
          exp(beta1[3] + cc*time.scale*beta1[nvar-1] + beta1[nvar] + gamma1*v) *
            integrate(lam1.fn, lower=tt, upper=tt+s)$value
        )) - surv.prob
      }
    }
    exit.time <- tt + get.sojourn.time(uniroot.fn)
    exit.idx <- as.integer(ceiling(exit.time * J)) # discretize exit time
    
    if (exit.idx > max.idx) {
      is.exit <- F
      exit.state <- NA
    } else {
      is.exit <- T
      
      # Determine exit state
      if (E==1) {
        x[2] <- (cc + exit.time - tt) * time.scale
      } else if (E==0 && Z==0) {
        r <- get.x4.idx(exit.time)
        x4 <- rep(0, RB)
        if (r <= RB) {
          x4[r] <- cc*time.scale
        }
        x[(1:RB)+3] <- x4
      }
      
      P <- P.fn(exit.time, x, v)
      diag(P) <- 0 # prevent "transitions" where they stay in the same state
      exit.state <- which(rmultinom(n=1, size=1, prob=P[state,]) > 0)
    }
  } else if (method == "discrete") {
    if (E == -1) {
      P.exit <- t(sapply(j:max.idx, function(u) {
        P.fn(u/J, x, v)[state,]
      }))
    } else if (E == 1) {
      P.exit <- t(sapply(j:max.idx, function(u) {
        x[2] <- cc+(u-j+1)/J
        P.fn(u/J, x, v)[state,]
      }))
    } else if (E == 0) {
      P.exit <- t(sapply(j:max.idx, function(u) {
        r <- findInterval(b.start + u - j/J, b.breaks)
        x4 <- rep(0, RB)
        if (r <= RB) {
          x4[r] <- cc*time.scale
        }
        x[(1:RB)+3] <- x4
        P.fn(u/J, x, v)[state,]
      }))
    }
    exit.states <- my.rmultinom(P.exit)
    
    if (any(exit.states != state)) {
      is.exit <- T
      exit.idx <- j-1 + min(which(exit.states != state)) # next time when not in current state
      exit.state <- exit.states[exit.idx - (j-1)]        # state occupied at that time
    } else {
      is.exit <- F
      exit.idx <- exit.state <- NA
    }
  }
  
  return(list(is.exit=is.exit, exit.idx=exit.idx, exit.state=exit.state))
}

generate.data <- function(method, n=nn, is.discrete.surv=F, print.increment=100) {
  print(paste0(Sys.time(), ": Generating v ..."))
  v.full <- rbinom(n, size=1, prob=pv)
  
  print(paste0(Sys.time(), ": Generating first exit times ..."))
  if (method == "cts") {
    surv.probs <- runif(n)
    if (is.discrete.surv) {
      expon <- 5
      div <- 10^expon
      surv.key.filename <- paste0("surv_key_e", expon, ".rds")
      if (file.exists(surv.key.filename)) {
        surv.key <- readRDS(surv.key.filename)
      } else {
        print(paste0(Sys.time(), ": Making survival probability key ..."))
        surv.key <- data.frame(
          v = rep(0:1, each=div+1),
          surv.idx = rep(0:div, times=2),
          surv.prob = rep((0:div)/div, times=2)
        )
        surv.key$first.exit.time <- get.first.exit.time(surv.key$v, surv.key$surv.prob)
        surv.key$first.exit.idx <- as.integer(ceiling(surv.key$first.exit.time * J))
        saveRDS(surv.key, file=surv.key.filename)
      }
      
      surv.df <- data.frame(v=v.full, surv.idx=as.integer(round(surv.probs, expon) * div))
      surv.df <- merge(surv.df, surv.key)
      first.exit.time <- surv.df$first.exit.time
      first.exit.idx <- surv.df$first.exit.idx
    } else {
      first.exit.time <- get.first.exit.time(v.full, surv.probs)
      first.exit.idx <- as.integer(ceiling(first.exit.time * J))
    }
  } else if (method == "discrete") {
    P.onset0 <- t(sapply(1:((A+tau)*J), function(u) { P.fn(u/J, x=rep(0,nvar), v=0)[1,] }))
    P.onset1 <- t(sapply(1:((A+tau)*J), function(u) { P.fn(u/J, x=rep(0,nvar), v=1)[1,] }))
    
    first.exit.idx <- rep(NA, n)
    first.exit.state <- rep(NA, n)
    
    sapply(1:n, function(i) {
      if (v.full[i]==0) {
        P <- P.onset0
      } else if (v.full[i]==1) {
        P <- P.onset1
      }
      onset.states <- my.rmultinom(P)
      
      if (any(onset.states != 1)) {
        first.exit.idx[i] <<- min(which(onset.states!=1)) # 1st time when not in (0^\circ, 0)
        first.exit.state[i] <<- onset.states[first.exit.idx[i]] # state occupied at that time
      } else {
        first.exit.idx[i] <<- (A+tau)*J + 1
      }
    })
    first.exit.time <- first.exit.idx / J
  }
  
  print(paste0(Sys.time(), ": Generating ages at selection time ..."))
  select.age <- runif(n, min=A-R.birth, max=A-L.birth)
  max.idx.full <- as.integer(ceiling((select.age + tau) * J)) # max number of observations for each i
  
  print(paste0(Sys.time(), ": Selecting the sample ..."))
  # Identify who exited (0^\circ, 0) by the selection time, by end of study, or neither
  i.exit <- which(first.exit.time <= select.age)
  i.late.exit <- setdiff(which(first.exit.time <= select.age + tau), i.exit)
  i.noexit <- setdiff(1:n, c(i.exit, i.late.exit))
  
  # # For those who did not exit (0^\circ, 0) by selection, randomly determine whether to select them
  # is.select0 <- rep(NA, n)
  # is.select0[-i.exit] <- as.logical(rbinom(n=n-length(i.exit), size=1, prob=p0))
  
  if (method == "cts") {
    # Identify which state was entered first after leaving (0^\circ, 0)
    first.exit.state <- rep(NA, n)
    first.exit.state <- get.first.exit.state(first.exit.time, v.full)
    # first.exit.state[i.exit] <- get.first.exit.state(first.exit.time[i.exit], v.full[i.exit])
    # 
    # i.late.exit.select <- intersect(i.late.exit, which(is.select0))
    # first.exit.state[i.late.exit.select] <- get.first.exit.state(
    #   first.exit.time[i.late.exit.select], v.full[i.late.exit.select]
    # )
  }
  
  # i.smoke <- which(first.exit.state == 5 & first.exit.time <= select.age)
  # i.select <- sort(c(which(is.select0), i.smoke))
  i.select <- 1:n
  
  print(paste0(Sys.time(), ": Generating follow-up data ..."))
  dat <- as.data.frame(data.table::rbindlist(lapply(i.select, function(i) {
    if (i==1 || i %% print.increment == 0) {
      print(paste0(Sys.time(), ": ", i))
    }
    
    max.idx <- max.idx.full[i]
    v <- v.full[i]
    dat.i <- data.frame(
      i = i,
      select.age = select.age[i],
      v = v,
      j = 1:max.idx,
      u = (1:max.idx)/J,
      E = NA,
      Z = NA,
      state.prev = NA,
      state = NA,
      c = NA,
      b = NA
    )
    
    fill.dat <- function(idx, E, Z, state.prev, state, cc, b) {
      dat.i$E[idx] <<- E
      dat.i$Z[idx] <<- Z
      dat.i$state.prev[idx] <<- state.prev
      dat.i$state[idx] <<- state
      dat.i$c[idx] <<- cc
      dat.i$b[idx] <<- b
    }
    
    if (i %in% i.noexit) {
      fill.dat(1:max.idx, E=-1, Z=0, state.prev=1, state=1, c=0, b=Inf)
    } else {
      fill.dat(1:(first.exit.idx[i]-1), E=-1, Z=0, state.prev=1, state=1, c=0, b=Inf)
      
      state <- first.exit.state[i]
      E <- get.E.for.state(state)
      Z <- get.Z.for.state(state)
      cc <- 0
      b <- ifelse(E==-1, Inf, 0)
      fill.dat(first.exit.idx[i], E, Z, state.prev=1, state=state, c=cc, b=b)
      
      j <- first.exit.idx[i] + 1
      
      # # For those who start smoking by selection, we still need to know if they're in (1,0) at selection
      # if (i %in% i.smoke) {
      #   select.j <- ceiling(select.age[i] * J)
      #   is.select <- F
      # } else {
      #   is.select <- T
      # }
      is.select <- T
      
      while (Z %in% c(0,10) && j <= max.idx) { # Z=10 means Z=1'
        exit.info <- get.exit.info(j, max.idx, cc, b, v, state, method)
        exit.idx <- exit.info$exit.idx
        
        if (!exit.info$is.exit) {
          break
        } else {
          exit.state <- exit.info$exit.state
          
          if (E == -1) { # never smoked
            if (exit.idx > j) {
              fill.dat(
                j:(exit.idx-1), E=E, Z=Z, state.prev=state, state=state, c=cc, b=b
              )
            }
          } else if (E == 1) { # currently smoking
            if (exit.idx > j) {
              fill.dat(
                j:(exit.idx-1),
                E=E, Z=Z, state.prev=state, state=state,
                c = cc + seq(1/J, by=1/J, length.out=exit.idx-j),
                b = 0
              )
            }
            cc <- cc + (exit.idx - j)/J
            b <- ifelse(get.E.for.state(exit.state)==1, b, 0)
          } else if (E == 0) { # previously smoked
            if (exit.idx > j) {
              fill.dat(
                j:(exit.idx-1),
                E=E, Z=Z, state.prev=state, state=state, c=cc,
                b = b + seq(1/J, by=1/J, length.out=exit.idx-j)
              )
            }
            b <- b + (exit.idx - j)/J
          }
          
          E <- get.E.for.state(exit.state)
          Z <- get.Z.for.state(exit.state)
          fill.dat(exit.idx, E, Z, state.prev=state, state=exit.state, c=cc, b=b)
          
          state <- exit.state
          j <- exit.idx + 1
        }
        
        # Check state occupied at selection time
        if (!is.select && j-1 >= select.j) {
          select.state <- dat.i$state[select.j]
          
          if (select.state == 5) { # (1,0)
            is.select <- ifelse(p1==1, T, as.logical(rbinom(n=1, size=1, prob=p1)))
          }
          if (!is.select) { return(NULL) }
        }
      }
      if (!is.select) { return(NULL) }
      
      # Current state is never exited; fill in remaining times accordingly
      if (j <= max.idx) {
        c.vec <- cc
        b.vec <- b
        if (Z %in% c(0,10)) { # unfailed
          if (E==1) {
            c.vec <- c.vec + seq(1/J, by=1/J, length.out=max.idx-j+1)
          } else if (E==0) {
            b.vec <- b.vec + seq(1/J, by=1/J, length.out=max.idx-j+1)
          }
        }
        fill.dat(j:max.idx, E, Z, state.prev=state, state=state, c=c.vec, b=b.vec)
      }
    }
    
    dat.i$E.prev <- c(-1, dat.i$E[1:(max.idx-1)])
    dat.i$Z.prev <- c(0, dat.i$Z[1:(max.idx-1)])
    
    return(dat.i)
  })))
  
  dat$B <- factor(findInterval(dat$b, b.breaks))
  dat$u.prev <- dat$u - 1/J
  dat$r <- dat$j - as.integer(ceiling(dat$select.age * J))
  
  dat$Z1 <- as.numeric(dat$Z==1)
  dat$Z2 <- as.numeric(dat$Z==2)
  
  dat$x2 <- ifelse(dat$E.prev == 1, dat$c*time.scale, 0)
  # dat.x4 <- as.data.frame(sapply(1:RB, function(r) {
  #   ifelse(dat$E.prev==0 & dat$Z.prev==0 & dat$B==r, dat$c*time.scale, 0)
  # }))
  # names(dat.x4) <- paste0("x4.", 1:RB)
  # dat <- cbind(dat, dat.x4)
  dat$x4 <- ifelse(dat$E.prev == 0 & dat$Z.prev == 0, dat$c*time.scale, 0)
  dat$x5 <- ifelse(dat$E.prev == 0 & dat$Z.prev == 10, dat$c*time.scale, 0)
  dat$x6 <- as.numeric(dat$Z.prev == 10)
  
  return(dat)
}
