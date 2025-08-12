rm(list=ls())
#### 0. Initialization ####
# Load packages
require(abind)
require(truncnorm)
require(extrafont)
loadfonts(device = "win")
# Set working directory to "Simulation - growth model"


# Load workspaces
load("Workspaces/res_mp3_ll_mm_pmm.RData")
load("Workspaces/res_mp3_lh_mm_pmm.RData")
load("Workspaces/res_mp3_hl_mm_pmm.RData")
load("Workspaces/res_mp3_hh_mm_pmm.RData")

load("Workspaces/res_mp3_ll_nm_pmm.RData")
load("Workspaces/res_mp3_lh_nm_pmm.RData")
load("Workspaces/res_mp3_hl_nm_pmm.RData")
load("Workspaces/res_mp3_hh_nm_pmm.RData")



# Source variable definitions
source("Functions/Functions_evaluate.R")
source("Functions/Functions.R")
source("Functions/VariableDefinitions.R")

# Find simulations with error (e.g. due to nonconvergence)
excl_ll_mm_pmm <- which(unlist(lapply(res_ll_mp3_mm, is.null)))
excl_lh_mm_pmm <- which(unlist(lapply(res_lh_mp3_mm, is.null)))
excl_hl_mm_pmm <- which(unlist(lapply(res_hl_mp3_mm, is.null)))
excl_hh_mm_pmm <- which(unlist(lapply(res_hh_mp3_mm, is.null)))

excl_ll_nm_pmm <- which(unlist(lapply(res_ll_mp3_nm, is.null)))
excl_lh_nm_pmm <- which(unlist(lapply(res_lh_mp3_nm, is.null)))
excl_hl_nm_pmm <- which(unlist(lapply(res_hl_mp3_nm, is.null)))
excl_hh_nm_pmm <- which(unlist(lapply(res_hh_mp3_nm, is.null)))

# Indices of first 2000 simulations without error
incl_ll_mm_pmm <- which(!(1:nSim%in% excl_ll_mm_pmm))[1:nSim_realized]
incl_lh_mm_pmm <- which(!(1:nSim%in% excl_lh_mm_pmm))[1:nSim_realized]
incl_hl_mm_pmm <- which(!(1:nSim%in% excl_hl_mm_pmm))[1:nSim_realized]
incl_hh_mm_pmm <- which(!(1:nSim%in% excl_hh_mm_pmm))[1:nSim_realized]

incl_ll_nm_pmm <- which(!(1:nSim%in% excl_ll_nm_pmm))[1:nSim_realized]
incl_lh_nm_pmm <- which(!(1:nSim%in% excl_lh_nm_pmm))[1:nSim_realized]
incl_hl_nm_pmm <- which(!(1:nSim%in% excl_hl_nm_pmm))[1:nSim_realized]
incl_hh_nm_pmm <- which(!(1:nSim%in% excl_hh_nm_pmm))[1:nSim_realized]

# Check for correct exclusion
range(incl_ll_mm_pmm)
range(incl_lh_mm_pmm)
range(incl_hl_mm_pmm)
range(incl_hh_mm_pmm)

range(incl_ll_nm_pmm)
range(incl_lh_nm_pmm)
range(incl_hl_nm_pmm)
range(incl_hh_nm_pmm)

# Select non-erroneous data for evaluation
ll_mm_pmm <- res_ll_mp3_mm[incl_ll_mm_pmm]
lh_mm_pmm <- res_lh_mp3_mm[incl_lh_mm_pmm]
hl_mm_pmm <- res_hl_mp3_mm[incl_hl_mm_pmm]
hh_mm_pmm <- res_hh_mp3_mm[incl_hh_mm_pmm]

ll_nm_pmm <- res_ll_mp3_nm[incl_ll_nm_pmm]
lh_nm_pmm <- res_lh_mp3_nm[incl_lh_nm_pmm]
hl_nm_pmm <- res_hl_mp3_nm[incl_hl_nm_pmm]
hh_nm_pmm <- res_hh_mp3_nm[incl_hh_nm_pmm]

#### 1. Fixed effects ####
# Extract fixed effects
fix_ll_mm_pmm <- simplify2array(lapply(ll_mm_pmm, "[[", 1))
fix_lh_mm_pmm <- simplify2array(lapply(lh_mm_pmm, "[[", 1))
fix_hl_mm_pmm <- simplify2array(lapply(hl_mm_pmm, "[[", 1))
fix_hh_mm_pmm <- simplify2array(lapply(hh_mm_pmm, "[[", 1))

fix_ll_nm_pmm <- simplify2array(lapply(ll_nm_pmm, "[[", 1))
fix_lh_nm_pmm <- simplify2array(lapply(lh_nm_pmm, "[[", 1))
fix_hl_nm_pmm <- simplify2array(lapply(hl_nm_pmm, "[[", 1))
fix_hh_nm_pmm <- simplify2array(lapply(hh_nm_pmm, "[[", 1))

dimnames(fix_ll_mm_pmm) <- dimnames(fix_ll_nm_pmm) <-
        dimnames(fix_lh_mm_pmm) <- dimnames(fix_lh_nm_pmm) <-
        dimnames(fix_hl_mm_pmm) <- dimnames(fix_hl_nm_pmm) <-
        dimnames(fix_hh_mm_pmm) <- dimnames(fix_hh_nm_pmm) <-
        list(c("b0", "b1"),
             c("Com w AR", "Com w/o AR", "RE-IMPUTE", "APPEND"),
             c("t3", "t4", "t5"),
             c(paste("sim", 1:nSim_realized)))

# Compute bias (estimated - true)
bias_fix_ll_mm_pmm <- apply(apply(fix_ll_mm_pmm, 2:4, function(x) x-c(b0,b1_lo)),1:3, mean)
bias_fix_lh_mm_pmm <- apply(apply(fix_lh_mm_pmm, 2:4, function(x) x-c(b0,b1_hi)),1:3, mean)
bias_fix_hl_mm_pmm <- apply(apply(fix_hl_mm_pmm, 2:4, function(x) x-c(b0,b1_lo)),1:3, mean)
bias_fix_hh_mm_pmm <- apply(apply(fix_hh_mm_pmm, 2:4, function(x) x-c(b0,b1_hi)),1:3, mean)

bias_fix_ll_nm_pmm <- apply(apply(fix_ll_nm_pmm, 2:4, function(x) x-c(b0,b1_lo)),1:3, mean)
bias_fix_lh_nm_pmm <- apply(apply(fix_lh_nm_pmm, 2:4, function(x) x-c(b0,b1_hi)),1:3, mean)
bias_fix_hl_nm_pmm <- apply(apply(fix_hl_nm_pmm, 2:4, function(x) x-c(b0,b1_lo)),1:3, mean)
bias_fix_hh_nm_pmm <- apply(apply(fix_hh_nm_pmm, 2:4, function(x) x-c(b0,b1_hi)),1:3, mean)

# Bias of third analysis (t=5)
bias_fix_ll_mm_pmm_t3 <- bias_fix_ll_mm_pmm[,,3]
bias_fix_lh_mm_pmm_t3 <- bias_fix_lh_mm_pmm[,,3]
bias_fix_hl_mm_pmm_t3 <- bias_fix_hl_mm_pmm[,,3]
bias_fix_hh_mm_pmm_t3 <- bias_fix_hh_mm_pmm[,,3]

bias_fix_ll_nm_pmm_t3 <- bias_fix_ll_nm_pmm[,,3]
bias_fix_lh_nm_pmm_t3 <- bias_fix_lh_nm_pmm[,,3]
bias_fix_hl_nm_pmm_t3 <- bias_fix_hl_nm_pmm[,,3]
bias_fix_hh_nm_pmm_t3 <- bias_fix_hh_nm_pmm[,,3]

# Bias of second analysis (t=4)
bias_fix_ll_mm_pmm_t2 <- bias_fix_ll_mm_pmm[,,2]
bias_fix_lh_mm_pmm_t2 <- bias_fix_lh_mm_pmm[,,2]
bias_fix_hl_mm_pmm_t2 <- bias_fix_hl_mm_pmm[,,2]
bias_fix_hh_mm_pmm_t2 <- bias_fix_hh_mm_pmm[,,2]

bias_fix_ll_nm_pmm_t2 <- bias_fix_ll_nm_pmm[,,2]
bias_fix_lh_nm_pmm_t2 <- bias_fix_lh_nm_pmm[,,2]
bias_fix_hl_nm_pmm_t2 <- bias_fix_hl_nm_pmm[,,2]
bias_fix_hh_nm_pmm_t2 <- bias_fix_hh_nm_pmm[,,2]

# Print bias
round(bias_fix_ll_mm_pmm,3)
round(bias_fix_lh_mm_pmm,3)
round(bias_fix_hl_mm_pmm,3)
round(bias_fix_hh_mm_pmm,3)

round(bias_fix_ll_nm_pmm,3)
round(bias_fix_lh_nm_pmm,3)
round(bias_fix_hl_nm_pmm,3)
round(bias_fix_hh_nm_pmm,3)

# Report range of bias
range.bias_fix <- range(bias_fix_ll_mm_pmm,
                        bias_fix_lh_mm_pmm,
                        bias_fix_hl_mm_pmm,
                        bias_fix_hh_mm_pmm,
                        bias_fix_ll_nm_pmm,
                        bias_fix_lh_nm_pmm,
                        bias_fix_hl_nm_pmm,
                        bias_fix_hh_nm_pmm)
#### 2. Random effects ####
# Extract random effects
ran_ll_mm_pmm <- simplify2array(lapply(ll_mm_pmm, "[[", 2))
ran_lh_mm_pmm <- simplify2array(lapply(lh_mm_pmm, "[[", 2))
ran_hl_mm_pmm <- simplify2array(lapply(hl_mm_pmm, "[[", 2))
ran_hh_mm_pmm <- simplify2array(lapply(hh_mm_pmm, "[[", 2))

ran_ll_nm_pmm <- simplify2array(lapply(ll_nm_pmm, "[[", 2))
ran_lh_nm_pmm <- simplify2array(lapply(lh_nm_pmm, "[[", 2))
ran_hl_nm_pmm <- simplify2array(lapply(hl_nm_pmm, "[[", 2))
ran_hh_nm_pmm <- simplify2array(lapply(hh_nm_pmm, "[[", 2))

dimnames(ran_ll_mm_pmm) <- dimnames(ran_ll_nm_pmm) <-
        dimnames(ran_lh_mm_pmm) <- dimnames(ran_lh_nm_pmm) <-
        dimnames(ran_hl_mm_pmm) <- dimnames(ran_hl_nm_pmm) <-
        dimnames(ran_hh_mm_pmm) <- dimnames(ran_hh_nm_pmm) <-
        list(c("U0", "U1", "e"),
             c("Com w AR", "Com w/o AR", "RE-IMPUTE", "APPEND"),
             c("t3", "t4", "t5"),
             c(paste("sim", 1:nSim_realized)))
# Compute bias random effects (estimated - true for U0, U1, e)
bias_ran_ll_mm_pmm <- round(apply(apply(ran_ll_mm_pmm, 2:4, function(x) x-c(U0,U1_lo,s2y)),1:3, mean),3)
bias_ran_lh_mm_pmm <- round(apply(apply(ran_lh_mm_pmm, 2:4, function(x) x-c(U0,U1_hi,s2y)),1:3, mean),3)
bias_ran_hl_mm_pmm <- round(apply(apply(ran_hl_mm_pmm, 2:4, function(x) x-c(U0,U1_lo,s2y)),1:3, mean),3)
bias_ran_hh_mm_pmm <- round(apply(apply(ran_hh_mm_pmm, 2:4, function(x) x-c(U0,U1_hi,s2y)),1:3, mean),3)

bias_ran_ll_nm_pmm <- round(apply(apply(ran_ll_nm_pmm, 2:4, function(x) x-c(U0,U1_lo,s2y)),1:3, mean),3)
bias_ran_lh_nm_pmm <- round(apply(apply(ran_lh_nm_pmm, 2:4, function(x) x-c(U0,U1_hi,s2y)),1:3, mean),3)
bias_ran_hl_nm_pmm <- round(apply(apply(ran_hl_nm_pmm, 2:4, function(x) x-c(U0,U1_lo,s2y)),1:3, mean),3)
bias_ran_hh_nm_pmm <- round(apply(apply(ran_hh_nm_pmm, 2:4, function(x) x-c(U0,U1_hi,s2y)),1:3, mean),3)
# Bias of third analysis (t=5)
bias_ran_ll_mm_pmm_t3 <- bias_ran_ll_mm_pmm[,,3]
bias_ran_lh_mm_pmm_t3 <- bias_ran_lh_mm_pmm[,,3]
bias_ran_hl_mm_pmm_t3 <- bias_ran_hl_mm_pmm[,,3]
bias_ran_hh_mm_pmm_t3 <- bias_ran_hh_mm_pmm[,,3]

bias_ran_ll_nm_pmm_t3 <- bias_ran_ll_nm_pmm[,,3]
bias_ran_lh_nm_pmm_t3 <- bias_ran_lh_nm_pmm[,,3]
bias_ran_hl_nm_pmm_t3 <- bias_ran_hl_nm_pmm[,,3]
bias_ran_hh_nm_pmm_t3 <- bias_ran_hh_nm_pmm[,,3]

# Bias of second analysis (t=4_)
bias_ran_ll_mm_pmm_t2 <- bias_ran_ll_mm_pmm[,,2]
bias_ran_lh_mm_pmm_t2 <- bias_ran_lh_mm_pmm[,,2]
bias_ran_hl_mm_pmm_t2 <- bias_ran_hl_mm_pmm[,,2]
bias_ran_hh_mm_pmm_t2 <- bias_ran_hh_mm_pmm[,,2]

bias_ran_ll_nm_pmm_t2 <- bias_ran_ll_nm_pmm[,,2]
bias_ran_lh_nm_pmm_t2 <- bias_ran_lh_nm_pmm[,,2]
bias_ran_hl_nm_pmm_t2 <- bias_ran_hl_nm_pmm[,,2]
bias_ran_hh_nm_pmm_t2 <- bias_ran_hh_nm_pmm[,,2]

# Print bias over t3-t5
round(bias_ran_ll_mm_pmm,3)
round(bias_ran_lh_mm_pmm,3)
round(bias_ran_hl_mm_pmm,3)
round(bias_ran_hh_mm_pmm,3)

round(bias_ran_ll_nm_pmm,3)
round(bias_ran_lh_nm_pmm,3)
round(bias_ran_hl_nm_pmm,3)
round(bias_ran_hh_nm_pmm,3)

# Report range of bias U1
range_bias_U1 <- range(bias_ran_ll_mm_pmm[2,3:4,],
                       bias_ran_lh_mm_pmm[2,3:4,],
                       bias_ran_hl_mm_pmm[2,3:4,],
                       bias_ran_hh_mm_pmm[2,3:4,],
                       bias_ran_ll_nm_pmm[2,3:4,],
                       bias_ran_lh_nm_pmm[2,3:4,],
                       bias_ran_hl_nm_pmm[2,3:4,],
                       bias_ran_hh_nm_pmm[2,3:4,])

# Report range of bias U0
range_bias_U0 <- range(bias_ran_ll_mm_pmm[1,3:4,],
                       bias_ran_lh_mm_pmm[1,3:4,],
                       bias_ran_hl_mm_pmm[1,3:4,],
                       bias_ran_hh_mm_pmm[1,3:4,],
                       bias_ran_ll_nm_pmm[1,3:4,],
                       bias_ran_lh_nm_pmm[1,3:4,],
                       bias_ran_hl_nm_pmm[1,3:4,],
                       bias_ran_hh_nm_pmm[1,3:4,])

# Report range of bias e
range_bias_et <- range(bias_ran_ll_mm_pmm[3,3:4,],
                       bias_ran_lh_mm_pmm[3,3:4,],
                       bias_ran_hl_mm_pmm[3,3:4,],
                       bias_ran_hh_mm_pmm[3,3:4,],
                       bias_ran_ll_nm_pmm[3,3:4,],
                       bias_ran_lh_nm_pmm[3,3:4,],
                       bias_ran_hl_nm_pmm[3,3:4,],
                       bias_ran_hh_nm_pmm[3,3:4,])
#### 3. Correlation ####
# Compute bias in correlation (imputed - complete)
bias_cor_ll_mm_pmm <- bias_cor(ll_mm_pmm)
bias_cor_lh_mm_pmm <- bias_cor(lh_mm_pmm)
bias_cor_hl_mm_pmm <- bias_cor(hl_mm_pmm)
bias_cor_hh_mm_pmm <- bias_cor(hh_mm_pmm)

bias_cor_ll_nm_pmm <- bias_cor(ll_nm_pmm)
bias_cor_lh_nm_pmm <- bias_cor(lh_nm_pmm)
bias_cor_hl_nm_pmm <- bias_cor(hl_nm_pmm)
bias_cor_hh_nm_pmm <- bias_cor(hh_nm_pmm)

# Bias third analysis (t=5)
bias_cor_ll_mm_pmm_t3 <- bias_cor_ll_mm_pmm[[3]]
bias_cor_lh_mm_pmm_t3 <- bias_cor_lh_mm_pmm[[3]]
bias_cor_hl_mm_pmm_t3 <- bias_cor_hl_mm_pmm[[3]]
bias_cor_hh_mm_pmm_t3 <- bias_cor_hh_mm_pmm[[3]]

bias_cor_ll_nm_pmm_t3 <- bias_cor_ll_nm_pmm[[3]]
bias_cor_lh_nm_pmm_t3 <- bias_cor_lh_nm_pmm[[3]]
bias_cor_hl_nm_pmm_t3 <- bias_cor_hl_nm_pmm[[3]]
bias_cor_hh_nm_pmm_t3 <- bias_cor_hh_nm_pmm[[3]]

# Bias second analysis (t=4)
bias_cor_ll_mm_pmm_t2 <- bias_cor_ll_mm_pmm[[2]]
bias_cor_lh_mm_pmm_t2 <- bias_cor_lh_mm_pmm[[2]]
bias_cor_hl_mm_pmm_t2 <- bias_cor_hl_mm_pmm[[2]]
bias_cor_hh_mm_pmm_t2 <- bias_cor_hh_mm_pmm[[2]]

bias_cor_ll_nm_pmm_t2 <- bias_cor_ll_nm_pmm[[2]]
bias_cor_lh_nm_pmm_t2 <- bias_cor_lh_nm_pmm[[2]]
bias_cor_hl_nm_pmm_t2 <- bias_cor_hl_nm_pmm[[2]]
bias_cor_hh_nm_pmm_t2 <- bias_cor_hh_nm_pmm[[2]]

# Bias first analysis (t=1)
bias_cor_ll_mm_pmm_t1 <- bias_cor_ll_mm_pmm[[1]]
bias_cor_lh_mm_pmm_t1 <- bias_cor_lh_mm_pmm[[1]]
bias_cor_hl_mm_pmm_t1 <- bias_cor_hl_mm_pmm[[1]]
bias_cor_hh_mm_pmm_t1 <- bias_cor_hh_mm_pmm[[1]]

bias_cor_ll_nm_pmm_t1 <- bias_cor_ll_nm_pmm[[1]]
bias_cor_lh_nm_pmm_t1 <- bias_cor_lh_nm_pmm[[1]]
bias_cor_hl_nm_pmm_t1 <- bias_cor_hl_nm_pmm[[1]]
bias_cor_hh_nm_pmm_t1 <- bias_cor_hh_nm_pmm[[1]]

# Bias in autocorrelation of outcome variables after imputation - first analysis (t=3)
Bias_autocor_ll_mm_t1 <- apply(apply(bias_cor_ll_mm_pmm_t1, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)]), 2, range)
Bias_autocor_lh_mm_t1 <- apply(apply(bias_cor_lh_mm_pmm_t1, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)]), 2, range)
Bias_autocor_hl_mm_t1 <- apply(apply(bias_cor_hl_mm_pmm_t1, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)]), 2, range)
Bias_autocor_hh_mm_t1 <- apply(apply(bias_cor_hh_mm_pmm_t1, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)]), 2, range)

Bias_autocor_ll_nm_t1 <- apply(apply(bias_cor_ll_nm_pmm_t1, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)]), 2, range)
Bias_autocor_lh_nm_t1 <- apply(apply(bias_cor_lh_nm_pmm_t1, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)]), 2, range)
Bias_autocor_hl_nm_t1 <- apply(apply(bias_cor_hl_nm_pmm_t1, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)]), 2, range)
Bias_autocor_hh_nm_t1 <- apply(apply(bias_cor_hh_nm_pmm_t1, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)]), 2, range)

# Bias in autocorrelation of outcome variables after imputation - second analysis (t=4)
Bias_autocor_ll_mm_t2 <- apply(apply(bias_cor_ll_mm_pmm_t2, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)]), 2, range)
Bias_autocor_lh_mm_t2 <- apply(apply(bias_cor_lh_mm_pmm_t2, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)]), 2, range)
Bias_autocor_hl_mm_t2 <- apply(apply(bias_cor_hl_mm_pmm_t2, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)]), 2, range)
Bias_autocor_hh_mm_t2 <- apply(apply(bias_cor_hh_mm_pmm_t2, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)]), 2, range)

Bias_autocor_ll_nm_t2 <- apply(apply(bias_cor_ll_nm_pmm_t2, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)]), 2, range)
Bias_autocor_lh_nm_t2 <- apply(apply(bias_cor_lh_nm_pmm_t2, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)]), 2, range)
Bias_autocor_hl_nm_t2 <- apply(apply(bias_cor_hl_nm_pmm_t2, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)]), 2, range)
Bias_autocor_hh_nm_t2 <- apply(apply(bias_cor_hh_nm_pmm_t2, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)]), 2, range)

# Bias in autocorrelation of outcome variables after imputation - third analysis (t=5)
Bias_autocor_ll_mm_t3 <- apply(apply(bias_cor_ll_mm_pmm_t3, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)]), 2, range)
Bias_autocor_lh_mm_t3 <- apply(apply(bias_cor_lh_mm_pmm_t3, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)]), 2, range)
Bias_autocor_hl_mm_t3 <- apply(apply(bias_cor_hl_mm_pmm_t3, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)]), 2, range)
Bias_autocor_hh_mm_t3 <- apply(apply(bias_cor_hh_mm_pmm_t3, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)]), 2, range)

Bias_autocor_ll_nm_t3 <- apply(bias_cor_ll_nm_pmm_t3, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)])
Bias_autocor_lh_nm_t3 <- apply(bias_cor_lh_nm_pmm_t3, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)])
Bias_autocor_hl_nm_t3 <- apply(bias_cor_hl_nm_pmm_t3, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)])
Bias_autocor_hh_nm_t3 <- apply(bias_cor_hh_nm_pmm_t3, 3, function(x) x[seq(nrow(x)+4,length(x), nrow(x)*2+2)])

# Print bias second analysis (t=4)
round(bias_cor_ll_mm_pmm_t2,3)
round(bias_cor_lh_mm_pmm_t2,3)
round(bias_cor_hl_mm_pmm_t2,3)
round(bias_cor_hh_mm_pmm_t2,3)

round(bias_cor_ll_nm_pmm_t2,3)
round(bias_cor_lh_nm_pmm_t2,3)
round(bias_cor_hl_nm_pmm_t2,3)
round(bias_cor_hh_nm_pmm_t2,3)

# Print bias third analysis (t=5)
round(bias_cor_ll_mm_pmm_t3,3)
round(bias_cor_lh_mm_pmm_t3,3)
round(bias_cor_hl_mm_pmm_t3,3)
round(bias_cor_hh_mm_pmm_t3,3)

round(bias_cor_ll_nm_pmm_t3,3)
round(bias_cor_lh_nm_pmm_t3,3)
round(bias_cor_hl_nm_pmm_t3,3)
round(bias_cor_hh_nm_pmm_t3,3)

# Report range for various combinations
range.bias_cor_mm_reimp_hi <- range(Bias_autocor_hl_mm_t3[,1],
                                 Bias_autocor_hh_mm_t3[,1])

range.bias_cor_nm_reimp_hi <- range(Bias_autocor_hl_nm_t3[,1],
                                 Bias_autocor_hh_nm_t3[,1])

range.bias_cor_mm_imp_hi <- range(Bias_autocor_hl_mm_t3[,2],
                               Bias_autocor_hh_mm_t3[,2])

range.bias_cor_nm_imp_hi <- range(Bias_autocor_hl_nm_t3[,2],
                               Bias_autocor_hh_nm_t3[,2])

range.bias_cor_mm_reimp_lo <- range(Bias_autocor_ll_mm_t3[,1],
                                    Bias_autocor_lh_mm_t3[,1])

range.bias_cor_nm_reimp_lo <- range(Bias_autocor_ll_nm_t3[,1],
                                    Bias_autocor_lh_nm_t3[,1])

range.bias_cor_mm_imp_lo <- range(Bias_autocor_ll_mm_t3[,2],
                                  Bias_autocor_lh_mm_t3[,2])

range.bias_cor_nm_imp_lo <- range(Bias_autocor_ll_nm_t3[,2],
                                  Bias_autocor_lh_nm_t3[,2])
#### 4. Coverage ####
# Extract coverage
cov_ll_mm_pmm <- coverage(ll_mm_pmm, b0, b1_lo)
cov_lh_mm_pmm <- coverage(lh_mm_pmm, b0, b1_hi)
cov_hl_mm_pmm <- coverage(hl_mm_pmm, b0, b1_lo)
cov_hh_mm_pmm <- coverage(hh_mm_pmm, b0, b1_hi)

cov_ll_nm_pmm <- coverage(ll_nm_pmm, b0, b1_lo)
cov_lh_nm_pmm <- coverage(lh_nm_pmm, b0, b1_hi)
cov_hl_nm_pmm <- coverage(hl_nm_pmm, b0, b1_lo)
cov_hh_nm_pmm <- coverage(hh_nm_pmm, b0, b1_hi)

# Extract coverage of slope (fixed effect)
cov_ll_mm_pmm_b1 <- t(cov_ll_mm_pmm[2,,])
cov_lh_mm_pmm_b1 <- t(cov_lh_mm_pmm[2,,])
cov_hl_mm_pmm_b1 <- t(cov_hl_mm_pmm[2,,])
cov_hh_mm_pmm_b1 <- t(cov_hh_mm_pmm[2,,])

cov_ll_nm_pmm_b1 <- t(cov_ll_nm_pmm[2,,])
cov_lh_nm_pmm_b1 <- t(cov_lh_nm_pmm[2,,])
cov_hl_nm_pmm_b1 <- t(cov_hl_nm_pmm[2,,])
cov_hh_nm_pmm_b1 <- t(cov_hh_nm_pmm[2,,])

# Print coverage slope
round(cov_ll_mm_pmm_b1,3)
round(cov_lh_mm_pmm_b1,3)
round(cov_hl_mm_pmm_b1,3)
round(cov_hh_mm_pmm_b1,3)

round(cov_ll_nm_pmm_b1,3)
round(cov_lh_nm_pmm_b1,3)
round(cov_hl_nm_pmm_b1,3)
round(cov_hh_nm_pmm_b1,3)

# Report range coverage
range.coverage <- range(cov_ll_mm_pmm_b1,
                        cov_lh_mm_pmm_b1,
                        cov_hl_mm_pmm_b1,
                        cov_hh_mm_pmm_b1,
                        cov_ll_nm_pmm_b1,
                        cov_lh_nm_pmm_b1,
                        cov_hl_nm_pmm_b1,
                        cov_hh_nm_pmm_b1)
#### 5. Efficiency ####
# Compute efficiency (of confidence interval width)
eff_ll_mm_pmm <- efficiency(ll_mm_pmm)
eff_lh_mm_pmm <- efficiency(lh_mm_pmm)
eff_hl_mm_pmm <- efficiency(hl_mm_pmm)
eff_hh_mm_pmm <- efficiency(hh_mm_pmm)

eff_ll_nm_pmm <- efficiency(ll_nm_pmm)
eff_lh_nm_pmm <- efficiency(lh_nm_pmm)
eff_hl_nm_pmm <- efficiency(hl_nm_pmm)
eff_hh_nm_pmm <- efficiency(hh_nm_pmm)

# Print efficiency
round(eff_ll_mm_pmm,3)
round(eff_lh_mm_pmm,3)
round(eff_hl_mm_pmm,3)
round(eff_hh_mm_pmm,3)

round(eff_ll_nm_pmm,3)
round(eff_lh_nm_pmm,3)
round(eff_hl_nm_pmm,3)
round(eff_hh_nm_pmm,3)

# Report efficiency
range.eff <- range(eff_ll_mm_pmm,
                   eff_lh_mm_pmm,
                   eff_hl_mm_pmm,
                   eff_hh_mm_pmm,
                   eff_ll_nm_pmm,
                   eff_lh_nm_pmm,
                   eff_hl_nm_pmm,
                   eff_hh_nm_pmm)
