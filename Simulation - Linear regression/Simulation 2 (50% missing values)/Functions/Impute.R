impute <-function(data, m, maxit = maxit){
  imp <- mice(data = data, m = m, maxit = maxit, print = FALSE, method = "norm")
  return(complete(imp, action = "long", include = TRUE))
}

#### Imputation, analysis and pooling reimputation ####
reimpute <- function(population.mu, population.lm, data.t1, data.t2, m, maxit){
  #no.covs <- ncol(data.t2)
  
  # Create storage objects
  #covered.mu <- covered.lm <- rep(0, ncol(data.t2))
  q.hat.mu <- u.mu <- fit <- #imp.sigma <- imp.st.sigma <- 
    vector("list", m)
  
  # Impute t2 data
  colnames(data.t2) <- c("x1", "y1", "x2", "y2")
  t2.imp <- mice(data = data.t2, m = m, maxit = maxit, print = FALSE, method = "norm")
  t.2 <- complete(t2.imp, action = "long", include = TRUE)
  
  # Analysis per dataset 
  for (dataset in 1:m){
    imp.data <- t.2[t.2$.imp == dataset,3:6]
    q.hat.mu[[dataset]] <- colMeans(imp.data)
    u.mu[[dataset]] <- cov(imp.data)/n
    fit[[dataset]] <- lm(formula = y2 ~ x1 + y1 + x2, data = imp.data)
    #imp.sigma[[dataset]] <- cov(imp.data[, 3:6])
    #imp.st.sigma[[dataset]] <- cor(imp.data[, 3:6])
    }
  
  # Pool results: Variable means
  q.bar.mu <- apply(simplify2array(q.hat.mu), 1, mean)
  var.u.mu <- diag(apply(simplify2array(u.mu), 1:2, mean))
  e.b.mu <- t(simplify2array(lapply(q.hat.mu, function(x) x - q.bar.mu)))
  b.mu <- t(e.b.mu) %*% e.b.mu * 1/(m-1)
  var.b.mu <- diag(b.mu)
  var.t.mu <- var.u.mu + (1 + 1/m) * var.b.mu
  df.mu <- df.function(var.u = var.u.mu, var.b = var.b.mu, var.t = var.t.mu, m1 = m, m2 = 1, n = n, no.covs = ncol(data.t2), no.parameters = 1)
  qt.mu <- qt(c(0.025, 0.975), df.mu)
  ci.95.mu <- cbind(q.bar.mu + qt.mu[1] * sqrt(var.t.mu), q.bar.mu + qt.mu[2] * sqrt(var.t.mu))

  # Pool results: Regression coefficients
  summary.lm <- as.data.frame(summary(pool(as.mira(fit)), conf.int=TRUE))
  q.bar.lm <- summary.lm$est
  ci.95.lm <- cbind(summary.lm[, "2.5 %"], summary.lm[, "97.5 %"])
  
  # Output
  output <- list("q.bar.mu"=q.bar.mu, "ci.95.mu"=ci.95.mu, 
                 "q.bar.lm"=q.bar.lm, "ci.95.lm"=ci.95.lm)
  return(output)
}








#### Imputation, analysis and pooling nested (NI) and appended (AI) imputation ####
nested.impute <- function(data.t1, population.mu, population.lm, new.data, m1, m2, maxit, n, pred){
 
  # Create storage objects
  no.covs <- ncol(data.t1) + ncol(new.data)
  #covered.mu <- covered.lm <- rep(0, ncol(mis.data))
  q.hat.mu <- q.hat.lm <- q.nest.mu <- q.nest.lm <- w.nest.mu <- w.nest.lm <- vector("list", m1)
  u.mu <- u.lm <- imp.st.sigma <-  vector("list", m1*m2)
  q.hat.nest.mu <- q.hat.nest.lm <- u.nest.mu <- u.nest.lm <- vector("list", m2)
  condition.m2 <- (m2>1)
  
   # Impute t1 data 
  t1.data <- impute(data.t1, m = m1, maxit = maxit)
  new.data.long <- new.data[rep(1:nrow(new.data), times = m1 + 1),]
  t2.data <- cbind(t1.data, new.data.long)

  # Impute t2-data, analyze and pool results
  for (imp.m1 in 1:m1){

    # Impute t2 data
    imp.t2 <- mice(data = t2.data[as.integer(t2.data$.imp) == imp.m1,], m = m2, maxit = maxit, pred = pred, print = FALSE, method = "norm")
    #if(!is.null(imp.t2$loggedEvents)){write.table(imp.t2$loggedEvents, file='mice_log.txt', append=TRUE) }
    
    # Analysis within nests
    for(imp.m2 in 1:m2){
      imp.data.nest <- complete(imp.t2, action = imp.m2)
      
      # Analysis: Variable means
      q.hat.nest.mu[[imp.m2]] <- colMeans(imp.data.nest[, 3:6])
      u.mu[[(imp.m1-1)*m2+imp.m2]] <- cov(imp.data.nest[, 3:6])/n
      
      # Analysis: Regression coefficients
      fit <- lm(formula = y2 ~ x1 + y1 + x2, data = imp.data.nest)
      q.hat.nest.lm[[imp.m2]] <- fit$coefficients
      u.lm[[(imp.m1-1)*m2+imp.m2]] <- vcov(fit)
      #imp.st.sigma[[(imp.m1-1)*m2+imp.m2]] <- cor(imp.data.nest[,3:6])
      }

      # Pooling within nests: Variable means
    q.hat.mu[[imp.m1]] <- q.hat.nest.mu
    #q.hat.nest.mu.mat <- simplify2array(q.hat.nest.mu)
    q.nest.mu[[imp.m1]] <- apply(simplify2array(q.hat.nest.mu), 1, mean)
    e.w.mu <- t(simplify2array(lapply(q.hat.nest.mu, function(x) x - q.nest.mu[[imp.m1]])))
    w.nest.mu[[imp.m1]] <- t(e.w.mu) %*% e.w.mu

    # Pooling within nests: Regression coefficients
    q.hat.lm[[imp.m1]] <- q.hat.nest.lm
    #q.hat.nest.lm.mat <- simplify2array(q.hat.nest.lm)
    q.nest.lm[[imp.m1]] <- apply(simplify2array(q.hat.nest.lm), 1, mean)
    e.w.lm <- t(simplify2array(lapply(q.hat.nest.lm, function(x) x - q.nest.lm[[imp.m1]])))
    w.nest.lm[[imp.m1]] <- t(e.w.lm) %*% e.w.lm
    }
  
  # Pooling between nests: Variable means
  q.bar.mu <- apply(simplify2array(q.nest.mu), 1, mean)
  var.u.mu <- diag(apply(simplify2array(u.mu), 1:2, mean))
  e.b.mu <- t(simplify2array(lapply(q.nest.mu, function(x) x - q.bar.mu)))
  var.b.mu <- diag(t(e.b.mu) %*% e.b.mu * m2/(m1-1))
  if(condition.m2) {var.w.mu <- 1 / (m1 * (m2 - 1)) * diag(apply(simplify2array(w.nest.mu), 1:2, sum))}
  else {var.w.mu <- rep(0, no.covs)}
  var.t.mu <- var.u.mu + 1/m2 * (1 + 1/m1) * var.b.mu + (1 - 1/m2) * var.w.mu

  df.mu <- df.function(var.u.mu, var.b.mu, var.w.mu, var.t.mu, m1, m2, n, no.covs, 1, c(3, 4))
  qt.mu <- qt(0.975, df.mu)
  ci.95.mu <- as.matrix(cbind(q.bar.mu - qt.mu * sqrt(var.t.mu), q.bar.mu + qt.mu * sqrt(var.t.mu)))

 
  
  # Pooling between nests: Regression coefficients
  q.bar.lm <- apply(simplify2array(q.nest.lm), 1, mean)
  var.u.lm <- diag(apply(simplify2array(u.lm), c(1:2), mean))
  e.lm <- t(simplify2array(lapply(q.nest.lm, function(x) x - q.bar.lm)))
  var.b.lm <- diag(t(e.lm)%*%e.lm * m2/(m1-1))
  if(condition.m2) var.w.lm <- 1/(m1*(m2-1)) * diag(apply(simplify2array(w.nest.lm), c(1:2), sum))
  else{var.w.lm <- rep(0, no.covs)}
  var.t.lm <- var.u.lm + 1/m2 * (1 + 1/m1) * var.b.lm + (1 - 1/m2) * var.w.lm
  
  df.lm <- df.function(var.u.lm, var.b.lm, var.w.lm, var.t.lm, m1, m2, n, no.covs, 4, c(3,4))
  qt.lm <- qt(0.975, df.lm)
  ci.95.lm <- as.matrix(cbind(q.bar.lm - qt.lm * sqrt(var.t.lm), q.bar.lm + qt.lm * sqrt(var.t.lm)))
 
  # Prepare output
  output <- list("q.bar.mu"=q.bar.mu, "ci.95.mu"=ci.95.mu, 
                 "q.bar.lm"=q.bar.lm, "ci.95.lm"=ci.95.lm)
                 #,"imp.st.sigma" = imp.st.sigma.mean)
  return(output)
}