rm(list=ls())

#### 0. Initalization ####
#### 0.1 set working directory ####
# Make sure folder "Simulation - Growth model" is set as working directory

#### 0.2 Load packages ####
require(mice) 
require(MASS)
require(mitml) 
require(miceadds) 
require(abind) 
require(MCMCglmm)
require(forecast)
require(truncnorm)
require(nlme) 

#### 0.2 Load functions and variable definitions####
source("Functions/Functions.R")
source("Functions/VariableDefinitions.R")
nSim <- nSim_realized <- 2
#### 1. Run simulation ####
#### 1.1 Monotone missingness pattern ####
set.seed(3000)
res_ll_mp3_mm <- simulate(J=J, cs=cs, rho_auto=rho_a1, rho_cross=rho_c1, nImp = nImp, mpat=mpat_mm, mperc=mperc, nSim=nSim, s2y=s2y, s2x=s2x, U0=U0,U1=U1_t,b0=b0, lag=lag, lag1=lag1, type="mm")
set.seed(3010)
res_lh_mp3_mm <- simulate(J=J, cs=cs, rho_auto=rho_a1, rho_cross=rho_c2, nImp = nImp, mpat=mpat_mm, mperc=mperc, nSim=nSim, s2y=s2y, s2x=s2x, U0=U0,U1=U1_t,b0=b0, lag=lag, lag1=lag1, type="mm")
set.seed(3020)
res_hl_mp3_mm <- simulate(J=J, cs=cs, rho_auto=rho_a2, rho_cross=rho_c1, nImp = nImp, mpat=mpat_mm, mperc=mperc, nSim=nSim, s2y=s2y, s2x=s2x, U0=U0,U1=U1_t,b0=b0, lag=lag, lag1=lag1, type="mm")
set.seed(3030)
res_hh_mp3_mm <- simulate(J=J, cs=cs, rho_auto=rho_a2, rho_cross=rho_c2, nImp = nImp, mpat=mpat_mm, mperc=mperc, nSim=nSim, s2y=s2y, s2x=s2x, U0=U0,U1=U1_t,b0=b0, lag=lag, lag1=lag1, type="mm")

#### 1.1.1 Save results ####
save(res_ll_mp3_mm,file="Workspaces/res_mp3_ll_mm_pmm.RData")
save(res_lh_mp3_mm,file="Workspaces/res_mp3_lh_mm_pmm.RData")
save(res_hl_mp3_mm,file="Workspaces/res_mp3_hl_mm_pmm.RData")
save(res_hh_mp3_mm,file="Workspaces/res_mp3_hh_mm_pmm.RData")

#### 1.2. Non-monotone missingness pattern ####
set.seed(3000)
res_ll_mp3_nm <- simulate(J=J,cs=cs, rho_auto=rho_a1, rho_cross=rho_c1, nImp = nImp, mpat=mpat_nm, mperc=mperc, nSim=nSim, s2y=s2y, s2x=s2x, U0=U0,U1=U1_t,b0=b0, lag=lag, lag1=lag1, type="nm")
set.seed(3010)
res_lh_mp3_nm <- simulate(J=J,cs=cs, rho_auto=rho_a1, rho_cross=rho_c2, nImp = nImp, mpat=mpat_nm, mperc=mperc, nSim=nSim, s2y=s2y, s2x=s2x, U0=U0,U1=U1_t,b0=b0, lag=lag, lag1=lag1, type="nm")
set.seed(3020)
res_hl_mp3_nm <- simulate(J=J,cs=cs, rho_auto=rho_a2, rho_cross=rho_c1, nImp = nImp, mpat=mpat_nm, mperc=mperc, nSim=nSim, s2y=s2y, s2x=s2x, U0=U0,U1=U1_t,b0=b0, lag=lag, lag1=lag1, type="nm")
set.seed(3030)
res_hh_mp3_nm <- simulate(J=J,cs=cs, rho_auto=rho_a2, rho_cross=rho_c2, nImp = nImp, mpat=mpat_nm, mperc=mperc, nSim=nSim, s2y=s2y, s2x=s2x, U0=U0,U1=U1_t,b0=b0, lag=lag, lag1=lag1, type="nm")

#### 1.2.1 Save results ####
save(res_ll_mp3_nm,file="Workspaces/res_mp3_ll_nm_pmm.RData")
save(res_lh_mp3_nm,file="Workspaces/res_mp3_lh_nm_pmm.RData")
save(res_hl_mp3_nm,file="Workspaces/res_mp3_hl_nm_pmm.RData")
save(res_hh_mp3_nm,file="Workspaces/res_mp3_hh_nm_pmm.RData")