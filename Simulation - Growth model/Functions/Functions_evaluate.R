#### Functions for evaluation ####
# Function to compute bias in the correlation matrices
#
# Arguments: 
# - res: list of results as obtained via simulate()
#
# Output: 
# - list of arrays with bias for analyses timepoints 1-3 (t = 3,4,5)
bias_cor <- function(res){
  # Extract correlation matrices
  cors <- lapply(res, "[[", 4)
 
  cors_com_t1 <- lapply(lapply(cors, function(x) x[[1]]), "[[", 1)
  cors_com_t2 <- lapply(lapply(cors, function(x) x[[1]]), "[[", 2)
  cors_com_t3 <- lapply(lapply(cors, function(x) x[[1]]), "[[", 3)
  cors_reimp <- lapply(cors, function(x) x[[3]])
  cors_imp <- lapply(cors, function(x) x[[4]])
  
  # For imputed data at analysis 1-3 (t=3,4,5)
  t1_reimp <- lapply(cors_reimp, function(x) x[[1]])
  t2_reimp <- lapply(cors_reimp, function(x) x[[2]])
  t3_reimp <- lapply(cors_reimp, function(x) x[[3]])
  t1_imp <- lapply(cors_imp, function(x) x[[1]])
  t2_imp <- lapply(cors_imp, function(x) x[[2]])
  t3_imp <- lapply(cors_imp, function(x) x[[3]])
  
  # Compute bias in correlation matrix (imputed - complete)
  bias_t1_reimp <- bias_t1_imp <- array(NA, dim=c(dim(cors_com_t1[[1]]),length(res)), dimnames=list(dimnames(cors_com_t1[[1]])[[1]],dimnames(cors_com_t1[[1]])[[2]],paste("sim", 1:nSim_realized)))
  bias_t2_reimp <- bias_t2_imp <-array(NA, dim=c(dim(cors_com_t2[[1]]),length(res)), dimnames=list(dimnames(cors_com_t2[[1]])[[1]],dimnames(cors_com_t2[[1]])[[2]],paste("sim", 1:nSim_realized)))
  bias_t3_reimp <- bias_t3_imp <- array(NA, dim=c(dim(cors_com_t3[[1]]),length(res)), dimnames=list(dimnames(cors_com_t3[[1]])[[1]],dimnames(cors_com_t3[[1]])[[2]],paste("sim", 1:nSim_realized)))

  for(sim in 1:length(res)){
    bias_t1_reimp[,,sim] <-  apply(simplify2array(lapply(t1_reimp[[sim]], `-`, cors_com_t1[[sim]])), 1:2, mean)
    bias_t2_reimp[,,sim] <-  apply(simplify2array(lapply(t2_reimp[[sim]], `-`, cors_com_t2[[sim]])), 1:2, mean)
    bias_t3_reimp[,,sim] <-  apply(simplify2array(lapply(t3_reimp[[sim]], `-`, cors_com_t3[[sim]])), 1:2, mean)
    bias_t1_imp[,,sim] <-  apply(simplify2array(lapply(t1_imp[[sim]], `-`, cors_com_t1[[sim]])), 1:2, mean)
    bias_t2_imp[,,sim] <-  apply(simplify2array(lapply(t2_imp[[sim]], `-`, cors_com_t2[[sim]])), 1:2, mean)
    bias_t3_imp[,,sim] <-  apply(simplify2array(lapply(t3_imp[[sim]], `-`, cors_com_t3[[sim]])), 1:2, mean)
  }
  
  bias_t1 <- abind(bias_t1_reimp, bias_t1_imp, along=4)
  bias_t1 <- aperm(bias_t1, c(1,2,4,3))
  
  bias_t2 <- abind(bias_t2_reimp, bias_t2_imp, along=4)
  bias_t2 <- aperm(bias_t2, c(1,2,4,3))
  
  
  bias_t3 <- abind(bias_t3_reimp, bias_t3_imp, along=4)
  bias_t3 <- aperm(bias_t3, c(1,2,4,3))
  
  dimnames(bias_t1)[[3]] <- dimnames(bias_t2)[[3]] <- dimnames(bias_t3)[[3]] <- c("RE-IMPUTE", "APPEND")
  
  bias.t1 <- round(apply(simplify2array(bias_t1), 1:3, mean), 3)
  bias.t2 <- round(apply(simplify2array(bias_t2), 1:3, mean), 3)
  bias.t3 <- round(apply(simplify2array(bias_t3), 1:3, mean), 3)
  
  return(list(bias.t1,bias.t2,bias.t3))
}

#### Coverage ####
# Function to compute coverage of the 95% confidence intervals of fixed effects
#
# Arguments: 
# - res: list of results as obtained via simulate()
#
# Output: 
# - arrays with coverage proportions for analyses timepoints 1-3 (t = 3,4,5)
coverage <- function(res, b0, b1){
  ci <- simplify2array(lapply(res, "[[", 3))
  dimnames(ci) <- list(c("b0", "b1"),
                       c("Com w AR", "Com w/o AR", "RE-IMPUTE", "APPEND"),
                       c("t3", "t4", "t5"),
                       c("lo", "hi"),
                       c(paste("sim", 1:length(res))))
  coverage_lo <- apply(ci[,,,1,], 2:4, function(x) c(b0, b1) > x)
  coverage_hi <- apply(ci[,,,2,], 2:4, function(x) c(b0, b1) < x)
  coverage <- apply(coverage_lo & coverage_hi, 1:3, mean)
  return(coverage)
}

#### Efficiency ####
# Function to compute relative efficiency of the 95% confidence intervals of fixed effects (RE-IMPUTE vs APPEND)
#
# Arguments: 
# - res: list of results as obtained via simulate()
#
# Output: 
# - array with relative efficiencies of analyses 1-3 (t = 3,4,5)

efficiency <- function(res){
  ci <- simplify2array(lapply(res, "[[", 3))
  dimnames(ci) <- list(c("b0", "b1"),
                       c("Com w AR", "Com w/o AR", "RE-IMPUTE", "APPEND"),
                       c("t3", "t4", "t5"),
                       c("lo", "hi"),
                       c(paste("sim", 1:length(res))))
  width <- apply(apply(ci, c(1,2,3,5), function(x) x[2]-x[1]), 1:3, mean)
  eff <- round(apply(width, c(1,3), function(x) x[4]/x[3]), 3)
  return(eff)
}