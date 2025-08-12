#### Variable definitions ####
J <- 425 # Number of subjects
cs <- 5 # Cluster size (number of timepoints)

mperc <- 0.50 # Missingness percentage (# of cases)

nImp <- 5 # m (number of imputations)
nSim <- 2100 # number of simulations performed (more than reported, such that occasional non-converging samples could be skipped)
nSim_realized <- 2000 # number of simulations to be reported
s2y <- s2x <- 1 # residual variance of outcome y_t and covariate x_t

lag1 <- 3 # first timepoint at which imputation/analysis occurs
lag <- 1 # number of timepoints between next imputation/analysis

#### Create missingness patterns ####
mpat <- expand.grid(c(0,1),c(0,1), c(0,1), c(0,1), c(0,1)) # possible combinations of missingness in outcome variable
pat_mm <- mpat[apply(mpat, 1, function(x) check.mp(x, mm=TRUE)),] # put restrictions on possible combinations for monotone missingness pattern
pat_nm <- mpat[apply(mpat, 1, function(x) check.mp(x, mm=FALSE)),] # put restrictions on possible combinations for nonmonotone missingness pattern

mpat_mm <- cbind(matrix(1, nrow=nrow(pat_mm), ncol=cs), pat_mm)                   # combine with covariates
mpat_mm <- cbind(1,  mpat_mm[, c(matrix(1:ncol(mpat_mm), nrow = 2, byrow = TRUE))]) # alternate columns to meet order of data

mpat_nm <- cbind(matrix(1, nrow=nrow(pat_nm), ncol=cs), pat_nm)                   # combine with covariates
mpat_nm <- cbind(1,  mpat_nm[, c(matrix(1:ncol(mpat_nm), nrow = 2, byrow = TRUE))]) # alternate columns to meet order of data


# Effect parameters 
U0 <- 0.20 # random intercept variance
U1_t <- 0.20 # random slope varance
rho_a1 <- rho_c1 <- 0.10 # low autocorrelation (rho_a1) and low crosscorrelation (rho_c1)
rho_a2 <- rho_c2 <- 0.50 # high autocorrelation (rho_a2) and high crosscorrelation (rho_c2)

b0 <- 1 # intercept

t.vars <- vector("list", cs) # list of time variables (0,...,T-1)
for(i in 1:(cs)) t.vars[[i]] <- c(paste("x.t.", i-1, sep = ""), paste("y.t.", i-1, sep = ""))

# adjusted means and variances of truncated normal distribution of crosscorrelation (high and low), to ensure a random slope variance of U1
mean_rho_tn_lo <- rho_c1 + 0.005 
var_rho_tn_lo <- U1_t + 0.005
mean_rho_tn_hi <- rho_c2 + 0.14
var_rho_tn_hi <- U1_t + 0.07

# true values for high and low cross-correlations
b1_lo <- etruncnorm(-1,1,mean_rho_tn_lo,sqrt(var_rho_tn_lo))
U1_lo <- vtruncnorm(-1,1,mean_rho_tn_lo,sqrt(var_rho_tn_lo))
b1_hi <- etruncnorm(-1,1,mean_rho_tn_hi,sqrt(var_rho_tn_hi))
U1_hi <- vtruncnorm(-1,1,mean_rho_tn_hi,sqrt(var_rho_tn_hi))