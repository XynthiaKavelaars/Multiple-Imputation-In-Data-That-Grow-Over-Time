resample <- function(mu, sigma, n, nSim, MisPerc, m1, m2, maxit){

  # Create storage objects
  ri <- ni <- ai <- rep(list(vector("list", nSim)), 2)
  qbar.mu.ri <-  qbar.lm.ri <-  
    qbar.mu.ni <-  qbar.lm.ni <-  
    qbar.mu.ai <-  qbar.lm.ai <-  matrix(NA, nrow = nSim, ncol = 4)
  
  # True (populations) means and regression coefficients
  population.mu <- mu
  cov.x <- sigma[c(1:3), c(1:3)]
  cov.xy <- sigma[4, c(1:3)]
  beta.hat <- solve(cov.x, cov.xy)
  intercept <- population.mu[4] - population.mu[c(1:3)] %*% beta.hat
  population.lm <- c(intercept, beta.hat)

  pb <- txtProgressBar(min = 0, max = nSim, style = 3)

  # Repeat simulation
  for(sample in 1:nSim){
    Sys.sleep(0.01)
    setTxtProgressBar(pb, sample)

    # Simulate data
    data <- simulate.data(mu, sigma, n)
    
    # Run study for each missingness pattern
    for (pat in 1:length(mis.pat)){
    # Make missing
    mis.data <- ampute(data, prop=MisPerc, mech="MCAR", patterns=mis.pat[[pat]], bycases=TRUE, run=TRUE)$amp

    data.t1 <- mis.data[,c(1,2)]
    data.t2 <- mis.data[,c(3:4)]
    #colnames(added.data) <- c("x2", "y2")
    
    # Apply multistage multiple imputation techniques, run complete data models and pool results
    ri[[pat]][[sample]] <- reimpute(population.mu, population.lm, data.t1, mis.data, m = m1, maxit = maxit)
    ni[[pat]][[sample]] <- nested.impute(data.t1, population.mu, population.lm, mis.data[,3:4], m1 = m1, m2 = m2, 
                                      maxit = maxit, n = n, pred = pred)
    ai[[pat]][[sample]] <- nested.impute(data.t1, population.mu, population.lm, mis.data[,3:4], m1 = m1, m2 = 1, 
                                      maxit = maxit, n = n, pred = pred)
  }}
  close(pb)
  
  output.list <- vector("list", length(mis.pat))
  for(pat in 1:length(mis.pat)){
  methods <- list(ri[[pat]], ni[[pat]], ai[[pat]])
  results <- c("q.bar.mu", "ci.95.mu", "q.bar.lm", "ci.95.lm")
  
  # Average bias & width, and sum coverage over samples
  q.bar.mu <- ci.mu <- mean.mu <- #sd.mu <- 
    bias.mu <- coverage.mu <- width.mu <-
    q.bar.lm <- ci.lm <- mean.lm <- #sd.lm <- 
    bias.lm <- coverage.lm <- width.lm <-
    #sigma.obs <- bias.sigma <- 
    vector("list", length(methods))
     for(method in 1:length(methods)){
      q.bar.mu[[method]] <- t(sapply(methods[[method]], "[[", "q.bar.mu", simplify = "array"))
      ci.mu[[method]] <- sapply(methods[[method]], "[[", "ci.95.mu", simplify = "array")
      mean.mu[[method]] <- colMeans(q.bar.mu[[method]])
      #sd.mu[[method]] <- apply(q.bar.mu[[method]], 2, sd)
      bias.mu[[method]] <- mean.mu[[method]]-population.mu
      coverage.mu[[method]] <- rowSums(apply(ci.mu[[method]], 3, function(x) population.mu>x[,1] & population.mu<x[,2])/nSim*100)
      width.mu[[method]] <- rowMeans(apply(ci.mu[[method]], 3, function(x) x[,2]-x[,1]))
      
      q.bar.lm[[method]] <- t(sapply(methods[[method]], "[[", "q.bar.lm", simplify = "array"))
      ci.lm[[method]] <- sapply(methods[[method]], "[[", "ci.95.lm", simplify = "array")
      mean.lm[[method]] <- colMeans(q.bar.lm[[method]])
      #sd.lm[[method]] <- apply(q.bar.lm[[method]], 2, sd)
      bias.lm[[method]] <- mean.lm[[method]]-population.lm
      coverage.lm[[method]] <- rowSums(apply(ci.lm[[method]], 3, function(x) population.lm>x[,1] & population.lm<x[,2])/nSim*100)
      width.lm[[method]] <- rowMeans(apply(ci.lm[[method]], 3,  function(x) x[,2]-x[,1]))
      
      #sigma.obs[[method]] <- sapply(methods[[method]], "[[", "imp.st.sigma", simplify = "array")
      #sigma.obs.mean <- apply(sigma.obs[[method]], 1:2, mean)
      #bias.sigma[[method]] <- sigma.obs.mean-sigma
    }
 output.list[[pat]] <- list("q.bar.mu"=q.bar.mu, "ci.mu"=ci.mu, "mean.mu"=mean.mu, #"sd.mu"=sd.mu, 
                            "bias.mu"=bias.mu, "coverage.mu"=coverage.mu, "width.mu"=width.mu,
                            "q.bar.lm"=q.bar.lm, "ci.lm"=ci.lm, "mean.lm"=mean.lm, #"sd.lm"=sd.lm,
                            "bias.lm"=bias.lm, "coverage.lm"=coverage.lm, "width.lm"=width.lm)
                            #,"sigma.obs"=sigma.obs, "bias.sigma"=bias.sigma)
 lapply(output.list[[pat]], function(x) names(x) <- c("ri", "ni", "ai"))}
   
  assign("output", output.list, envir = .GlobalEnv)
  return(output.list)}

