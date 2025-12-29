
logit <- Vectorize(function(p) { log(p/(1-p)) })

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
lam20.fn <- Vectorize(function(t) {
  kap2*rho2 * (rho2*t)^(kap2-1)
  # same as: dweibull(t, kap2, 1/rho2) / pweibull(t, kap2, 1/rho2, lower.tail=F)
})

###### INCORPORATING SYMPTOMS #########

phi <- alpha1
delta <- beta1

omega <- logit(c(0.5, 0.75))

gamma0 <- log(8)
gamma1 <- log(4)
gamma2 <- c(0, log(2))

###### SIMULATION PARAMETERS ##########

R <- 1000  # size of partition of [0,tau)
n <- 1000 # sample size

nsim <- 1000 # number of simulation replicates




