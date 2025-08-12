#### 0. Simulation setup ####
# Make sure folder "Simulation 1 (20% Missingness)" is set as working directory
### Basics #####
require(MASS)#install.packages("MASS")
require(pwr)#install.packages("pwr")
require(mice)#install.packages("mice")
require(Matrix)
require(matrixcalc)#install.packages("matrixcalc")
require(doParallel)#install.packages("doParallel")

no.cores <- detectCores()
registerDoParallel(cores=no.cores-1)

#### Required sample size ####
test.data <- mvrnorm(n=100, rep(0,4), matrix(c(1, 0.1, 0.1, 0.1, 
                                               0.1, 1, 0.1, 0.1,
                                               0.1, 0.1, 1, 0.1,
                                               0.1, 0.1, 0.1, 1),
                                             nrow=4, ncol=4, byrow=TRUE), empirical=TRUE)
multiple.r.squared <- summary(lm(test.data[,4]~test.data[,1]+test.data[,2]+test.data[,3]))$r.squared
f.squared <- multiple.r.squared/(1-multiple.r.squared)

# Required sample size for 80 or 95% power respectively, under smallest effect size
n.80 <- pwr.f2.test(u = 3, v = NULL, f2 = f.squared, sig.level = 0.05, power = 0.80)$v 

#### Required number of samples: Pilot ####

mu <- rep(0, 4)
set.seed(1545)
n <- round(n.80)
n.sample <- 50
mis.perc <- 0.2
m1.1 <- 5
m1.2 <- 5
m2 <- 5
maxit <- 25

source("Functions/Condition.R")

rho.vector.sel <- c(0.1, 0.3, 0.5, 0.7)
save.vector <- condition.list <- rep(NA, length(rho.vector.sel)*length(rho.vector.sel))
for(rho.t1 in rho.vector.sel){
  for(rho.t2 in rho.vector.sel){
    filename <- paste0("out_", rho.t1*10, rho.t2*10)
    condition.name <- paste0("Conditions/r_0", rho.t1*100, "_0", rho.t2*100, ".R")
    cond <- (which(rho.vector.sel==rho.t1)-1)*length(rho.vector.sel)+which(rho.vector.sel==rho.t2)
    save.vector[cond] <- filename
    condition.list[cond] <- condition.name
  }
}

# Create conditions
condition.function <- function(condition.list, index, save.vector){
  source(condition.list[index])
  output <- get(save.vector[index])
  return(output)
}

# Run simulation
results <- foreach(cond=1:length(condition.list), 
                   .packages=c("mice", "MASS", "matrixcalc", "Matrix")) %dopar% {
                     # Load functions
                     source("Functions/Simulate.R")
                     source("Functions/Impute.R")
                     source("Functions/Resample.R")
                     source("Functions/Sigma.R")
                     source("Functions/Df.R")
                     condition.function(condition.list, cond, save.vector)
                   }
names(results) <- save.vector
list2env(results, envir=.GlobalEnv)
save(list=save.vector, file = "Workspaces/Simulation_design.RData")


#### Find sds ####
load("Workspaces/Simulation_design.RData")
rho.vector <- c(0.1, 0.3, 0.5, 0.7)
results.function <- function(type, rho.vector, mp){
  results <- array(NA, dim=c(length(rho.vector), length(rho.vector),4,3))
  i <- 1
  for(rho.t1 in 1:7){
    j <- 1
    condition.rho.t1 <- (rho.t1/10) %in% rho.vector 
    if(condition.rho.t1){
      for(rho.t2 in 1:7){
        condition.rho.t2 <- (rho.t2/10) %in% rho.vector
        if(condition.rho.t2){
          filename <- paste0("out_", rho.t1, rho.t2)
          for(var in 1:4){
            for(method in 1:3){
              results[i, j, var, method] <- get(filename)[[mp]][[type]][[method]][var]
            }}
          j <- j+1}
      }
      i <- i+1}}
  return(results)}

sds.mu.monotone <- results.function(type= "sd.mu", rho.vector=rho.vector, mp=1)
sds.mu.nonmonotone <- results.function(type= "sd.mu", rho.vector=rho.vector, mp=2)
sds.lm.monotone <- results.function(type= "sd.lm", rho.vector=rho.vector, mp=1)
sds.lm.nonmonotone <- results.function(type= "sd.lm", rho.vector=rho.vector, mp=2)

max.sds <- max(sds.mu.monotone, sds.mu.nonmonotone, sds.lm.monotone, sds.lm.nonmonotone)
mean.sds.mu.monotone <- mean(sds.mu.monotone)
mean.sds.mu.nonmonotone <- mean(sds.mu.nonmonotone)
mean.sds.lm.monotone <- mean(sds.lm.monotone)
mean.sds.lm.nonmonotone <- mean(sds.lm.nonmonotone)

mean.sd <- mean(mean.sds.mu.monotone, mean.sds.mu.nonmonotone, 
                mean.sds.lm.monotone, mean.sds.lm.nonmonotone)

par(mfrow=c(2,2))
hist(as.vector(sds.mu.monotone))
hist(as.vector(sds.mu.nonmonotone))
hist(as.vector(sds.lm.monotone))
hist(as.vector(sds.lm.nonmonotone))

# Max sd = 0.13; round up to 0.20
selected.sd <- 0.20

#### Required number of samples ####
# 5% significance, error of 0.01
n.sample <- (qnorm(0.975)*selected.sd/0.01)^2
