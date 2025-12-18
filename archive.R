
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
  
  P.onset00 <- t(sapply(1:R, function(r) { P.fn(u[r], cc=0, v=0, w=0)[1,] }))
  P.onset01 <- t(sapply(1:R, function(r) { P.fn(u[r], cc=0, v=0, w=1)[1,] }))
  P.onset10 <- t(sapply(1:R, function(r) { P.fn(u[r], cc=0, v=1, w=0)[1,] }))
  P.onset11 <- t(sapply(1:R, function(r) { P.fn(u[r], cc=0, v=1, w=1)[1,] }))
  
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
    
    # If the person ever leaves state (0^\circ, 0) ...
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
      c.prev <- 0
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


