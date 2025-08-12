### 1. Execute simulation #####
# Make sure folder "Simulation" is set as working directory
require(MASS)
require(mice)
require(Matrix)
require(matrixcalc)
require(doParallel)

no.cores <- detectCores()
registerDoParallel(cores=no.cores-1)

# Create conditions
source("Functions/GetVariableDefinitions.R")
source("Functions/Condition.R")


# Presented conditions
Rho <- c(0.1, 0.3, 0.5, 0.7)

#### Simulation ####
save.vector <- condition.list <- rep(NA, length(Rho)*length(Rho))
for(rho.t1 in Rho){
  for(rho.t2 in Rho){
    filename <- paste0("out_", rho.t1*10, rho.t2*10)
    condition.name <- paste0("Conditions/r_0", rho.t1*100, "_0", rho.t2*100, ".R")
    cond <- (which(Rho==rho.t1)-1)*length(Rho)+which(Rho==rho.t2)
    save.vector[cond] <- filename
    condition.list[cond] <- condition.name
  }
}

condition.function <- function(condition.list, index, save.vector){
  source(condition.list[index])
  output <- get(save.vector[index])
  return(output)
}

results <- foreach(cond=1:length(condition.list), 
                   .packages=c("mice", "MASS", "matrixcalc", "Matrix")) %dopar% {
                     # Load functions
                     source("Functions/Simulate.R")
                     source("Functions/Impute.R")
                     source("Functions/Resample.R")
                     source("Functions/Sigma.R")
                     source("Functions/Df.R")
                     source("Functions/GetVariableDefinitions.R")
                     condition.function(condition.list, cond, save.vector)
                       }
names(results) <- save.vector
list2env(results, envir=.GlobalEnv)
save(list=save.vector, file = "Workspaces/Simulation.RData")


