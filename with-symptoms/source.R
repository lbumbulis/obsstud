
logit <- Vectorize(function(p) { log(p/(1-p)) })

get.coverage <- Vectorize(function(est, se, target, conf.level=0.95) {
  crit <- qnorm(conf.level + (1-conf.level)/2)
  return(as.numeric(est - crit*se <= target && target <= est + crit*se))
}, vectorize.args=c("est", "se"))

#######################################
# time/age scale: 1 unit time = 100 years
time.scale <- 100

tau <- 50/time.scale # at most 50 years of follow-up from selection date
pv <- 1/2

RB <- 10 # number of intervals for quit backward recurrence time b(t)
b.breaks <- seq(0, (RB-1)/time.scale, length.out=RB) # breaks defining the intervals

nvar <- RB + 5 # number of parameters in each of alpha and beta

###### DISEASE PROCESS ################
## Covariate coefficients
beta1p <- log(c(1.2, 1.05, 1.05, 1.05^((RB-(1:RB))/RB), 1, 1))
gamma1p <- log(1.1)

beta1 <- log(c(1.2, 1.05, 1.05, 1.05^((RB-(1:RB))/RB), 1.05, 10)) # consider increasing the 10 to 20
gamma1 <- log(1.1)

beta2 <- log(c(1.05, 1.02, 1.02, rep(1,RB), 1.02, 1))
gamma2 <- log(1.2)

# Temporary simplification # TODO: remove
beta1p[2] <- beta1[2] <- beta2[2] <- 0
beta1p[4:(nvar-2)] <- beta1[4:(nvar-2)] <- 0
beta1[nvar-1] <- beta2[nvar-1] <- 0

## Baseline intensities
# Lung cancer-free death
lam2.multiroot.fn <- function(logkr) {
  k <- exp(logkr[1])
  r <- exp(logkr[2])
  t.med <- (log(2))^(1/k) / r
  t.90 <- (log(10))^(1/k) / r
  return(c(t.med - 75/time.scale, t.90 - 90/time.scale))
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
  F60 <- (1-pv)*integrate(integrand0, 0, 60/time.scale)$value +
    pv*integrate(integrand1, 0, 60/time.scale)$value
  F80 <- (1-pv)*integrate(integrand0, 0, 80/time.scale)$value +
    pv*integrate(integrand1, 0, 80/time.scale)$value
  
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

alpha1 <- log(c(1, 1, 1, rep(1,RB), 1, 1))
eta1 <- log(1.2)

alpha0 <- log(c(1, 1, 1, rep(1,RB), 1, 100)) # 10, 20 are too low
eta0 <- log(0.8)

## Baseline intensities
# Starting to smoke (1st time)
psi1c.multiroot.fn <- function(logkr) {
  k <- exp(logkr[1])
  r <- exp(logkr[2])
  
  psi1c.fn <- function(tt) { k*r * (r*tt)^(k-1) } # hazard
  Psi1c.fn <- function(tt) { (r*tt)^k }           # cumulative hazard
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
  FE10 <- (1-pv)*integrate(integrand0, 0, 10/time.scale)$value +
    pv*integrate(integrand1, 0, 10/time.scale)$value
  FE40 <- (1-pv)*integrate(integrand0, 0, 40/time.scale)$value +
    pv*integrate(integrand1, 0, 40/time.scale)$value
  
  return(c(FE10 - 0.01, FE40 - 1/3))
}
temp <- rootSolve::multiroot(psi1c.multiroot.fn, start=c(0, 0))
psi11 <- exp(temp$root[1]) # like kappa (shape)
psi12 <- exp(temp$root[2]) # like rho (1/scale)
psi1c.fn <- Vectorize(function(tt) {
  psi11*psi12 * (psi12*tt)^(psi11-1)
  # same as: dweibull(tt, psi11, 1/psi12) / pweibull(tt, psi11, 1/psi12, lower.tail=F)
})

# Quitting and resuming smoking
rho.E <- 3
psi0 <- 0.1 / (1+rho.E)
psi1 <- rho.E * psi0

###### SAMPLING MECHANISM #############
L.birth <- 0             # left endpoint of birth cohort (can think of as 1920)
R.birth <- 30/time.scale # right endpoint (1950)
A <- 70/time.scale       # sampling date (1990)

p0 <- 0.5 # prob of selecting individual who is in (0^\circ, 0) at selection time
p1 <- 1     # prob of selecting individual who is in (1,0) at selection time

###### SIMULATION PARAMETERS ##########
J <- 1000 # 1/J is the increment used to discretize time

nn <- 10^4   # super-sample size
nsim <- 1000 # number of simulation replicates


