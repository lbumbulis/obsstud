
args <- commandArgs(trailingOnly=TRUE)
start.iter <- as.numeric(args[1])
iter.jump <- as.numeric(args[2])
datagen.method <- args[3]
model.type <- args[4]

stop.iter <- start.iter + iter.jump - 1

if (model.type == "cox") {
  library(survival)
}

source("source.R")
source("datagen.R")

sim.seeds <- readRDS("sim_seeds_nsim1000.rds")
.Random.seed <- sim.seeds[[(start.iter %/% iter.jump) + 1]]
# For iter.jump=10: start.iter=1 -> 1st seed, 11 -> 2nd, 21 -> 3rd, ...



for (iter in start.iter:stop.iter) {
  print(paste0(Sys.time(), ": Starting iter ", iter))
  system.time(dat <- generate.data(datagen.method, print.increment=1))
  
  # ## Features of the data, to be saved later
  # v.full <- tapply(dat$v, dat$i, function(x) x[1])
  # 
  # # rows = time points, cols = states; used for plots
  # state.occupancy.j <- with(dat, table(j, state)) # time 0 is birth
  # state.occupancy.r <- with(dat, table(r, state)) # time 0 is recruitment
  # 
  # transition.count <- with(
  #   dat, table(state.prev, state, i)
  # ) # used for number of quit attempts per person, among smokers
  # 
  # save(
  #   v.full, state.occupancy.j, state.occupancy.r, transition.count,
  #   file = paste0("./data_features/data_features_iter", iter, ".RData")
  # )
  
  # if (iter==1) {
  #   saveRDS(dat, file="test_dat_iter1.rds", compress="bzip2")
  # }
  
  
  
  # Filter the data to what is needed for analysis
  dat.sub <- dat[-which(dat$state.prev==dat$state & dat$Z %in% 1:2),] # already in absorbing state
  dat.sub <- dat.sub[-which(dat.sub$Z.prev==10),] # not at risk for -> 2 once symptoms develop
  
  if (model.type == "cox") {
    fail.times <- unique(dat.sub$j[which(dat.sub$Z==2)])
    dat.cox <- dat.sub[which(dat.sub$j %in% fail.times),] # only need data from times when someone fails
    # # Check how many ties there are in the failure times
    # length(fail.times)
    # length(dat$j[which(dat$Z==2 & dat$Z.prev==0)])
    # hist(dat$j[which(dat$Z==2 & dat$Z.prev==0)], breaks=seq(0.4, 1.1, by=0.0001)*J, main="", xlab="j")
    
    print(paste0(Sys.time(), ": Starting Cox analysis"))
    m.cox <- coxph(
      Surv(u.prev, u, Z2) ~ factor(E.prev) + v,
      data = dat.cox, method = "breslow"
    )
    # method="breslow": takes < 10 seconds for J=1000 (n=10^4), < 1min for J=10^4 (n=5000)
    # method="exact": takes over 1.5h for J=1000 (n=10^4), over 5h for J=10^4 (n=5000)
    
    coef.est <- coef(m.cox)[c(2,1,3)]
    coef.var <- diag(vcov(m.cox))[c(2,1,3)]
    
    coef.est.filename <- "mcox_cause2_est.csv"
    coef.var.filename <- "mcox_cause2_var.csv"
    
    H0.hat <- basehaz(m.cox, centered=F)
    saveRDS(H0.hat, file=paste0("./basehaz/basehaz_cause2_iter", iter, ".rds"), compress="bzip2")
    
    # saveRDS(m.cox, file=paste0("./models/mcox_cause2_iter", iter, ".rds"), compress="bzip2")
  } else {
    dat.discrete <- aggregate(
      Z2 ~ j + E.prev + v, data = dat.sub,
      FUN = function(x) { c(y = sum(x), m = length(x)) }
    )
    dat.discrete <- do.call(data.frame, dat.discrete)
    
    if (model.type == "logis") {
      m.logis <- glm(
        cbind(Z2.y, Z2.m - Z2.y) ~ -1 + factor(j) + factor(E.prev) + v,
        data = dat.discrete, family = binomial
      )
      ncov <- length(coef(m.logis))
      coef.est <- coef(m.logis)[c(ncov-c(1,2,0), 1:(ncov-3))]
      coef.var <- diag(vcov(m.logis))[c(ncov-c(1,2,0), 1:(ncov-3))]
      
      coef.est.filename <- "mlogis_cause2_est.csv"
      coef.var.filename <- "mlogis_cause2_var.csv"
      
    } else if (model.type == "pois") {
      # m.pois <- glm(
      #   Z2.y ~ -1 + factor(j) + factor(E.prev) + v,
      #   data = dat.discrete, # also tried dat.sub
      #   family = poisson
      # )
    }
  }
  
  write.table(
    cbind(iter=iter, t(coef.est)), file=coef.est.filename,
    append=file.exists(coef.est.filename), quote=F, sep=",",
    row.names=F, col.names=!file.exists(coef.est.filename)
  )
  
  write.table(
    cbind(iter=iter, t(coef.var)), file=coef.var.filename,
    append=file.exists(coef.var.filename), quote=F, sep=",",
    row.names=F, col.names=!file.exists(coef.var.filename)
  )
}








