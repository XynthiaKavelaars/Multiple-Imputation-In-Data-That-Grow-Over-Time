#### Functions ####

#### 1. Check missingness pattern ####
# Function to check whether the missingnesspattern meets specific criteria, namely:
# - monotone or nonmonotone
# - at least 2 out of 3 first timepoints observed to ensure an identifiable imputation model
# - at least 1 variable unobserved

# Arguments:
# - x: vector of length ncol(data) with ones and zeros, representing a missingness pattern (0=unobserved; 1=observed)
# - mm: monotone (TRUE) or nonmonotone (FALSE)

# Output:
# check: Logical; if TRUE, missingness pattern meets criteria; if FALSE missingness pattern does not meet criteria
check.mp <- function(x, mm = TRUE){
  monotone <- !is.unsorted(rev(x)) # if monotone
  nonmonotone <- is.unsorted(rev(x)) # if nonmonotone 
  n_all_obs <- sum(x) < cs # not all variables observed
  f_s_obs <- x[1] & (x[2] | x[3]) # 2 out of 3 first timepoints observed
  if(mm) {check <- monotone & n_all_obs & f_s_obs
  }else{check <- nonmonotone & n_all_obs & f_s_obs}
  return(check)
}

#### 2. Ampute data ####
# Function to create missing data via MCAR-mechanism
#
#Arguments:
# - data: n x p matrix with complete data of n subjects and p variables (wide format)
# - prop: proportion of cases with incomplete data
# - patterns: m x p matrix with m missingness patterns of p variables

# Output:
# amp: n x p matrix with incomplete data

ampute_MCAR <- function(data, prop=0.5, patterns){
  n <- nrow(data)
  p <- ncol(data)
  npat <- nrow(patterns)
  
  nmis <- round(prop * n)
  mp <- sample(1:npat, nmis, replace=TRUE)
  idMis <- sample(1:n, nmis, replace=FALSE)
  
  amp <- data
  
  for(i in 1:nmis){
    dat.i <- data[idMis[i],]
    dat.i[patterns[mp[i],] == 0] <- NA
    amp[idMis[i],] <- dat.i
  }
  return(amp)
} 

#### 3. Create lag-1 autocorrelation matrix ####
# Function to create a lag-1 autocorrelation matrix
#
# Arguments: 
# - cs: number of timepoints
# - rho: lag 1 autocorrelation
#
# Output:
# - mat: cs x cs symmetric autocorrelation matrix 

ar1_cor <- function(cs, rho) {
  exponent <- abs(matrix(1:cs - 1, nrow = cs, ncol = cs, byrow = TRUE) - 
                    (1:cs - 1))
  mat <- rho^exponent
  return(mat)
}


#### 4. Simulate ####
# Function to perform simulation: generate, ampute, impute, and analyze data
#
# Arguments:
# - J: no. of subjects
# - cs: cluster size (number of timepoints T)
# - rho_auto: lag-1 autocorrelation
# - rho_cross: cross correlation between x_t and y_t
# - nImp: number of imputation for RE-IMPUTE and first lag of APPEND
# - mpat: m x p matrix with m missingness patterns of p variables
# - mperc: missingness proportion
# - nSim: number of simulation
# - s2y: residual variance in outcome variable y_t
# - s2x: residual variance in covariate x_t
# - U0: random intercept variance
# - U1: random slope variance
# - b0: fixed intercept effect
# - lag1: first timepoints of imputation and analysis
# - lag: number of timepoints between imputations/analyses
# - type: type of missingness pattern (mm = monotone; nm = nonmonotone)
# 
# Output:
# - res: list of length nSim with a nested list of estimated fixed effects, random effects, confidence intervals, correlation matrices and (pooled) estimates for complete and imputed data
simulate <- function(J, cs, rho_auto, rho_cross, nImp, mpat, mperc, nSim, s2y, s2x, U0, U1, b0, lag, lag1, type="mm"){
  # Create output file for error management
  outputFile <- paste0("mp_",mperc*10,"_ra_",rho_auto*10,"_rc_",rho_cross*10,"_", type, "_pmm.txt")
  
    ts <- rep(1:cs - 1, J) # Vector of timepoints
  
  pid <- seq_len(J)  # individual id
  pid <- rep(pid, each = cs)  # repeat each ID cs times
  
  # Adjust mean and variance of truncated normal distribution to ensure a random slope variance of U1
  if(rho_cross==rho_c1){
    mean_rho_tn <- rho_c1 + 0.005
    var_rho_tn <- U1 + 0.005
  } else if(rho_cross==rho_c2){
    mean_rho_tn <- rho_c2 + 0.14
    var_rho_tn <- U1 + 0.07
  }

# True values given truncated normal crosscorrelation
b1 <- etruncnorm(-1,1,mean_rho_tn,sqrt(var_rho_tn)) # Expected value
U1 <- vtruncnorm(-1,1,mean_rho_tn,sqrt(var_rho_tn)) # Variance

gam <- c(b0, b1) # vector of true fixed effects
G <- matrix(c(U0, 0, 0, U1), nrow=2) # matrix of random variances

Sigma_Y <- ar1_cor(cs, rho_auto) # autocorrelation matrix for outcome
Sigma_X <- ar1_cor(cs, rho_auto) # autocorrelation matrix for covariate

t.vars <- vector("list", cs) # list of time variables
for(i in 1:(cs)) t.vars[[i]] <- c(paste("x.t.", i-1, sep = ""), paste("y.t.", i-1, sep = ""))

t_analysis <- seq(lag1, cs, lag) # vector of timepoints of analyses

res <- vector("list", nSim)
for(sim in 1:nSim){
  skip_to_next <- FALSE # For later use wot tryCatch(). If TRUE: in case of error, go to next simulation
  x.t <- mvrnorm(J, rep(0, cs), Sigma_X, empirical=TRUE) # Sample covariate data x_1
  X_i <- cbind(1, seq_len(cs) - 1)  # for each individual
  X <- cbind(X_i[rep(1:nrow(X_i),times=J),], c(t(x.t)))  # repeat each row cs times
 
  uj <- mvrnorm(J, mu = rep(0, length(gam)), Sigma = G, empirical=TRUE) # Generate person-level (lv-2) random effects
  eij <- mvrnorm(J, mu = rep(0, cs), Sigma = Sigma_Y, empirical=TRUE) # Generate repeated-measure-level (lv-1) error term

  # Compute individual regression coefficients
  betaj <- matrix(gam, nrow = J, ncol = 2, byrow = TRUE) + uj
  betaj <- betaj[rep(1:nrow(betaj), each=cs),]
  
  # Compute outcome:
  y <- rowSums(X[,-2] * betaj)  +  c(t(eij))

  # Form a data frame
  dat.long <- cbind.data.frame(pid,ts, c(t(x.t)),y)
  colnames(dat.long) <- c("pid", "time", "x.t", "y.t")
  dat.wide <- reshape(dat.long, v.names=c("x.t", "y.t"), idvar="pid", timevar="time", direction="wide")

  # Ampute data
  amp.wide <- ampute_MCAR(data=dat.wide, prop=mperc, patterns=mpat)
  amp.long <- reshape(amp.wide, idvar="pid", varying=unlist(t.vars), v.names=c("x.t", "y.t"), times=c(0:(cs-1)), direction="long")
  amp.long.sort <- amp.long[order(amp.long[,1], amp.long[,2]),]

  # Create predictor matrix for imputation
  pred <- make.predictorMatrix(amp.long.sort)
  pred <- matrix(0, nrow=ncol(amp.long.sort), ncol=ncol(amp.long.sort))
  colnames(pred) <- rownames(pred) <- colnames(amp.long.sort)
  pred["y.t", ] <- c(-2, 0, 2, 2)
  
  # Create storage objects
  reimp_pmm <- reimplist_pmm <- 
    imp_pmm <- implist_pmm <-  
     vector("list", length(t_analysis))
  res.com_AR <- res.com <-  
    dat_pmm_wide <- dat_pmm_long <- reimp_pmm.wide <- imp_pmm.wide <- imp_pmm.long <-
    vector("list", length(t_analysis))
  
  # Impute first timepoint of analysis via RE-IMPUTE and APPEND
  reimp_pmm[[1]] <- imp_pmm[[1]] <- tryCatch({
  mice(amp.long.sort[amp.long.sort$time %in% 0:(t_analysis[1]-1),], pred = pred, meth = "2l.pmm", m=nImp, print = FALSE, maxit = 20)
  },error = function(e) {
    write.table(paste0("Error in sim ", sim, " (re)imputation lag 1. Error ", as.character(e)) , outputFile, append=TRUE)
    skip_to_next <<- TRUE
    return(list())})
  if(skip_to_next) { next }  
  reimplist_pmm[[1]] <- implist_pmm[[1]] <- as.mitml.list(complete(imp_pmm[[1]], action="all", include=FALSE))  
  
  # Reshape long-wide
  vars <- unlist(t.vars[1:t_analysis[1]])
  dat_pmm_wide[[1]] <- amp.wide[,c("pid", unlist(t.vars[1:t_analysis[1]]))]
  dat_pmm_long[[1]] <- reshape(dat_pmm_wide[[1]], idvar="pid", varying=vars, v.names=c("x.t", "y.t"), times=0:(t_analysis[1]-1), direction="long")
  reimp_pmm.wide[[1]] <- imp_pmm.wide[[1]] <- lapply(implist_pmm[[1]], function(x) reshape(x, idvar="pid", v.names=c("x.t", "y.t"), timevar="time", direction="wide"))

  # Perform complete data analysis of first timepoint of imputation/analysis with explicit autocorrelation structure
    res.com_AR[[1]] <- tryCatch({
  lme(y.t ~ 1 + x.t, random = ~ 1 + x.t | pid, data = dat.long[dat.long$time %in% 0:(t_analysis[1]-1),], correlation = corAR1(form = ~ time), control=lmeControl(opt="optim", maxIter = 1e8, msMaxIter = 1e8))
  }, error = function(e) {
    write.table(paste0("Error in sim ", sim, " complete data analysis, autocorrelation, lag 1. Error ", as.character(e)) , outputFile, append=TRUE)
    skip_to_next <<- TRUE
    return(list())})
  if(skip_to_next) { next }   
 
    # Perform complete data analysis of first timepoint of imputation/analysis without explicit autocorrelation structure
  res.com[[1]] <- tryCatch({
    lme(y.t ~ 1 + x.t, random = ~ 1 + x.t | pid, data = dat.long[dat.long$time %in% 0:(t_analysis[1]-1),], control=lmeControl(opt="optim", maxIter = 1e8, msMaxIter = 1e8))
  }, error = function(e) {
    write.table(paste0("Error in sim ", sim, " complete data analysis, no autocorrelation, lag 1. Error ", as.character(e)) , outputFile, append=TRUE)
    skip_to_next <<- TRUE
    return(list())})
  if(skip_to_next) { next }
 
  # Perform imputation and analysis for all other timepoints
  for(i in t_analysis[-1]){
  ind <- which(t_analysis ==i)
  vars <- unlist(t.vars[1:i])
  new_vars <- c((t_analysis[ind-1]+1):i)
   dat_pmm_wide[[ind]] <- lapply(imp_pmm.wide[[ind-1]], function(x) cbind(x,amp.wide[,unlist(t.vars[i])]))
  dat_pmm_long[[ind]] <- lapply(dat_pmm_wide[[ind]], function(x) reshape(x, idvar="pid", varying=vars,v.names=c("x.t", "y.t"), times=c(0:(i-1)), direction="long"))

  # RE-IMPUTE
    reimp_pmm[[ind]] <- tryCatch({
    mice(amp.long.sort[amp.long.sort$time %in% 0:(i-1),], pred = pred, meth = "2l.pmm", m=nImp, print = FALSE, maxit = 20)
  },error = function(e) {
    write.table(paste0("Error in sim ", sim, ", reimputation, lag ", i, ". Error ", as.character(e)) , outputFile, append=TRUE)
    skip_to_next <<- TRUE
    return(list())})
  if(skip_to_next) { break }   
  reimplist_pmm[[ind]] <- as.mitml.list(complete(reimp_pmm[[ind]], action="all", include=FALSE))

    reimp_pmm.wide[[ind]] <- lapply(reimplist_pmm[[ind]], function(x) reshape(x, idvar="pid", v.names=c("x.t", "y.t"), timevar="time", direction="wide"))

  # APPEND
    # Check for incomplete data in new variables. If present, impute via APPEND 
   if(any(is.na(amp.wide[,unlist(t.vars[i])]))){ 
    imp_pmm[[ind]] <- tryCatch({
      lapply(dat_pmm_long[[ind]], function(x) mice(x, pred = pred, meth = "2l.pmm", m=1, print = FALSE, maxit = 20))
       },error = function(e) {
      write.table(paste0("Error in sim ", sim, "appended imputation, lag", i, ". Error ", as.character(e)) , outputFile, append=TRUE)
      skip_to_next <<- TRUE
      return(list())})
     if(skip_to_next) { break }   
       implist_pmm[[ind]] <- lapply(imp_pmm[[ind]], function(x) complete(x, action="all", include=FALSE)[[1]])
  
  imp_pmm.wide[[ind]] <- lapply(implist_pmm[[ind]], function(x) reshape(x, idvar="pid", v.names=c("x.t", "y.t"), timevar="time", direction="wide"))

  #If absent, concatenate new data without imputation
} else {
     imp_pmm.wide[[ind]] <- dat_pmm_wide[[ind]]
    implist_pmm[[ind]] <- dat_pmm_long[[ind]]

}
    
    # Perform analysis for all analysis timepoints on complete data - explicit autocorrelation structure
  res.com_AR[[ind]] <- tryCatch({
  lme(y.t ~ 1 + x.t, random = ~ 1 + x.t | pid, data = dat.long[dat.long$time %in% 0:(i-1),], correlation = corAR1(form = ~ time), control=lmeControl(opt="optim", maxIter = 1e8, msMaxIter = 1e8))
  }, error = function(e) {
    write.table(paste0("Error in sim ", sim, " complete data analysis, autocorrelation, lag", i, ". Error ", as.character(e)) , outputFile, append=TRUE)
    skip_to_next <<- TRUE
    return(list())})
  if(skip_to_next) { break }   
 
  # Perform analysis for all analysis timepoints on complete data - no explicit autocorrelation structure
   res.com[[ind]] <- tryCatch({
    lme(y.t ~ 1 + x.t, random = ~ 1 + x.t | pid, data = dat.long[dat.long$time %in% 0:(i-1),], control=lmeControl(opt="optim", maxIter = 1e8, msMaxIter = 1e8))
  }, error = function(e) {
    write.table(paste0("Error in sim ", sim, " complete data analysis, no autocorrelation, lag", i, ". Error ", as.character(e)) , outputFile, append=TRUE)
    skip_to_next <<- TRUE
    return(list())})
  if(skip_to_next) { break } 
  }
  if(skip_to_next) { next }   
  
  implist_pmm <- lapply(implist_pmm, as.mitml.list)
  reimplist_pmm <- lapply(reimplist_pmm, as.mitml.list)

  # Perform analysis on imputed data
  # RE-IMPUTE
  res.reimp_pmm <- tryCatch({
      lapply(reimplist_pmm, function(x) with(x, lme(y.t ~ 1 + x.t, random = ~ 1 + x.t | pid, correlation = corAR1(form = ~ time), control=lmeControl(opt="optim", maxIter = 1e8, msMaxIter = 1e8))))
  }, error = function(e) {
      write.table(paste0("Error in sim ", sim, " analysis after reimputation. Error ", as.character(e)) , outputFile, append=TRUE)
      skip_to_next <<- TRUE
      return(list())})
  if(skip_to_next) { next }
  
  # APPEND
  res.imp_pmm <- tryCatch({
  lapply(implist_pmm, function(x) with(x, lme(y.t ~ 1 + x.t, random = ~ 1 + x.t | pid, correlation = corAR1(form = ~ time), control=lmeControl(opt="optim", maxIter = 1e8, msMaxIter = 1e8))))
      }, error = function(e) {
        write.table(paste0("at sim ", sim, " analysis after appended imputation. Error ", as.character(e)) , outputFile, append=TRUE)
        skip_to_next <<- TRUE
        return(list())})
      if(skip_to_next) { next }  
  
  # Compute correlation matrices of imputed data
Cor.reimp_pmm <- lapply(reimp_pmm.wide, function(x) lapply(x, function(y) cor(y[,-1])))
Cor.imp_pmm <- lapply(imp_pmm.wide, function(x) lapply(x, function(y) cor(y[,-1])))

  # Collect estimates of fitted models
  Estimates.imp_pmm <- lapply(res.imp_pmm, testEstimates, var.comp=TRUE)
  Estimates.reimp_pmm <- lapply(res.reimp_pmm, testEstimates, var.comp=TRUE)
   estimates <- list(Estimates.reimp_pmm, Estimates.imp_pmm)
  
   # Obtain fixed effects
  Fix.com_AR <- simplify2array(lapply(res.com_AR, fixef))
  Fix.com <- simplify2array(lapply(res.com, fixef))
  Fix.reimp_pmm <- simplify2array(lapply(res.reimp_pmm, function(x) testEstimates(x)$estimates[,1]))
  Fix.imp_pmm <- simplify2array(lapply(res.imp_pmm, function(x) testEstimates(x)$estimates[,1]))
  fixefs <- abind(Fix.com_AR, Fix.com, Fix.reimp_pmm, Fix.imp_pmm, along=3)#, Fix.reimp_pmm, Fix.imp_pmm, Fix.reimp_norm, Fix.imp_norm, along=3)
  fixefs <- aperm(fixefs, c(1,3,2))
  
  # Obtain confidence intervals of fixed effects
  CI_Fix.com_AR <- simplify2array(lapply(res.com_AR, function(x) intervals(x, which="fixed")$fixed[,c(1,3)]))
  CI_Fix.com <- simplify2array(lapply(res.com, function(x) intervals(x, which="fixed")$fixed[,c(1,3)]))
  CI_Fix.imp_pmm <- simplify2array(lapply(Estimates.imp_pmm, function(x) confint(x, oldnames=FALSE)))
  CI_Fix.reimp_pmm <- simplify2array(lapply(Estimates.reimp_pmm, function(x) confint(x, oldnames=FALSE)))
  ci.fix <- abind(CI_Fix.com_AR, CI_Fix.com, CI_Fix.reimp_pmm, CI_Fix.imp_pmm, along=4)#, CI_Fix.reimp_pmm, CI_Fix.imp_pmm, CI_Fix.reimp_norm, CI_Fix.imp_norm, along=4)
  ci.fix <- aperm(ci.fix, c(1,4,3,2))

  # Obtain random effects of imputed data
  Ran.imp_pmm <- simplify2array(lapply(Estimates.imp_pmm, function(x) x$var.comp[c(1,3,4)]))
  Ran.reimp_pmm <- simplify2array(lapply(Estimates.reimp_pmm, function(x) x$var.comp[c(1,3,4)]))
 
  # Obtain correlation matrices and random effects of complete data
    Cor.com_AR <- Cor.com <- vector("list", length(t_analysis))
   Ran.com_AR <- Ran.com <- matrix(NA, 3, length(t_analysis))
  for(i in 1:(length(t_analysis))){
    Cor.com_AR[[i]] <- Cor.com[[i]] <- cor(dat.wide[,c(2:(t_analysis[i]*2+1))])
    Ran.com_AR[,i] <- as.numeric(VarCorr(res.com_AR[[i]])[,1])
   Ran.com[,i] <- as.numeric(VarCorr(res.com[[i]])[,1])}
  
  ranefs <- abind(Ran.com_AR, Ran.com, Ran.reimp_pmm, Ran.imp_pmm, along=3)#, Ran.reimp_pmm, Ran.imp_pmm, Ran.reimp_norm, Ran.imp_norm, along=3)
  ranefs <- aperm(ranefs, c(1,3,2))
  
  cors <- list(Cor.com_AR, Cor.com, Cor.reimp_pmm, Cor.imp_pmm)

   res[[sim]] <- 
     list(fixefs, ranefs, ci.fix, cors, estimates)#autocors)
}
return(res)
}

