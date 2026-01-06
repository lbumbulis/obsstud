
logit <- Vectorize(function(p) { log(p/(1-p)) })

tau <- 1
pv <- 1/2

RB <- 10       # number of intervals for quit backward recurrence time b(t)
b.breaks <- c(seq(0, ((RB-1)/10)/RB, length.out=RB), Inf) # breaks defining the intervals

nvar <- RB + 4 # number of parameters in each of alpha and beta

###### DISEASE PROCESS ################
## Covariate coefficients
beta1p <- log(c(1.2, 1.05, 1.05, 1+0.05*(RB-(1:RB)+1)/RB, 1))
gamma1p <- log(1.1)

beta1 <- log(c(1.2, 1.05, 1.05, 1+0.05*(RB-(1:RB)+1)/RB, 5)) # consider increasing the 5 to 10 or 20
gamma1 <- log(1.1)

beta2 <- log(c(1.05, 1.02, 1.02, rep(1,RB), 1))
gamma2 <- log(1.2)

## Baseline intensities
# Lung cancer-free death
lam2.multiroot.fn <- function(logkr) {
  k <- exp(logkr[1])
  r <- exp(logkr[2])
  t.med <- (log(2))^(1/k) / r
  t.90 <- (log(10))^(1/k) / r
  return(c(t.med - 0.75, t.90 - 0.9))
}
temp <- rootSolve::multiroot(lam2.multiroot.fn, start=c(0,0))
kap2 <- exp(temp$root[1])
rho2 <- exp(temp$root[2])
lam2.fn <- Vectorize(function(tt) {
  kap2*rho2 * (rho2*tt)^(kap2-1)
  # same as: dweibull(tt, kap2, 1/rho2) / pweibull(tt, kap2, 1/rho2, lower.tail=F)
})

# Lung cancer
rho <- 5 # multiplier from baseline lung cancer intensity to baseline symptom intensity
# consider increasing from 5 to 10

lam1.multiroot.fn <- function(logkr) {
  k <- exp(logkr[1])
  r <- exp(logkr[2])
  
  l1.fn <- function(tt) { k*r * (r*tt)^(k-1) }
  L1.fn <- function(tt) { (r*tt)^k } # same as -pweibull(t, k, 1/r, lower=F, log=T)
  Lam2.fn <- function(tt) { (rho2*tt)^kap2 }
  
  integrand0 <- Vectorize(function(s) {
    l1star <- l1.fn(s) * (1+rho)
    L1star <- L1.fn(s) * (1+rho)
    l1star * exp(-(L1star + Lam2.fn(s)))
  })
  integrand1 <- Vectorize(function(s) {
    l1star <- l1.fn(s) * (exp(gamma1) + rho*exp(gamma1p))
    L1star <- L1.fn(s) * (exp(gamma1) + rho*exp(gamma1p))
    l1star * exp(-(L1star + Lam2.fn(s)*exp(gamma2)))
  })
  F60 <- (1-pv)*integrate(integrand0, 0, 0.6)$value + pv*integrate(integrand1, 0, 0.6)$value
  F80 <- (1-pv)*integrate(integrand0, 0, 0.8)$value + pv*integrate(integrand1, 0, 0.8)$value
  
  return(c(F60 - 0.02, F80 - 0.05))
}
temp <- rootSolve::multiroot(lam1.multiroot.fn, start=c(1,-1))
kap1 <- exp(temp$root[1])
rho1 <- exp(temp$root[2])
lam1.fn <- Vectorize(function(tt) {
  kap1*rho1 * (rho1*tt)^(kap1-1)
  # same as: dweibull(tt, kap1, 1/rho1) / pweibull(tt, kap1, 1/rho1, lower.tail=F)
})

# Lung cancer symptoms
lam1p.fn <- Vectorize(function(tt) { rho*lam1.fn(tt) })

###### SMOKING PROCESS ################
## Covariate coefficients
eta1c <- log(1.2) # TODO: Check this

alpha1 <- log(c(1, 1, 1, rep(1,RB), 1))
eta1 <- log(1.2)

alpha0 <- log(c(1, 1, 1, rep(1,RB), 10)) # consider increasing the 10 to 100
eta0 <- log(0.8)

## Baseline intensities
# Starting to smoke (1st time)
psi1c.multiroot.fn <- function(logkr) {
  k <- exp(logkr[1])
  r <- exp(logkr[2])
  
  psi1c.fn <- function(tt) { k*r * (r*tt)^(k-1) }
  Psi1c.fn <- function(tt) { (r*tt)^k }
  Lam1.fn <- function(tt) { (rho1*tt)^kap1 }
  Lam2.fn <- function(tt) { (rho2*tt)^kap2 }
  
  integrand0 <- Vectorize(function(s) {
    Lam1star <- Lam1.fn(s) * (1+rho)
    psi1c.fn(s) * exp(-(Lam1star + Lam2.fn(s) + Psi1c.fn(s)))
  })
  integrand1 <- Vectorize(function(s) {
    Lam1star <- Lam1.fn(s) * (exp(gamma1) + rho*exp(gamma1p))
    psi1c.fn(s)*exp(eta1c) * exp(-(Lam1star + Lam2.fn(s)*exp(gamma2) + Psi1c.fn(s)*exp(eta1c)))
  })
  FE25 <- (1-pv)*integrate(integrand0, 0, 0.25)$value + pv*integrate(integrand1, 0, 0.25)$value
  FE70 <- (1-pv)*integrate(integrand0, 0, 0.7)$value + pv*integrate(integrand1, 0, 0.7)$value
  
  return(c(FE25 - 0.3, FE70 - 0.35))
}
temp <- rootSolve::multiroot(psi1c.multiroot.fn, start=c(0, 0))
psi11 <- exp(temp$root[1]) # like kappa
psi12 <- exp(temp$root[2]) # like rho
psi1c.fn <- Vectorize(function(tt) {
  psi11*psi12 * (psi12*tt)^(psi11-1)
  # same as: dweibull(tt, psi11, 1/psi12) / pweibull(tt, psi11, 1/psi12, lower.tail=F)
})

# Quitting and resuming smoking
rho.E <- 3
psi0 <- 0.1 / (1+rho.E)
psi1 <- rho.E * psi0

###### SIMULATION PARAMETERS ##########
R <- 1000 # size of partition of [0,tau)

nn <- 1000    # sample size
nsim <- 1000 # number of simulation replicates


