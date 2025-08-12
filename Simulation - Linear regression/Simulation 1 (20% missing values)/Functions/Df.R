
# Calculate degrees of freedom according to pooling rules Rubin (1987) and Rubin (2003)
df.function <- function(var.u, var.b, var.w = rep(0, no.covs), var.t, m1, m2 = 1, n, no.covs, 
                        no.parameters, t2.parameters = NULL) {
  
  df.old <- lambda <- df <- df.observed <- rep(NA, no.covs)
  condition.m2 <- (m2 > 1)
  
  for(variable in 1:no.covs){
    b.term <- (1/(m1-1)) * (1/m2 * (1 + 1/m1) * var.b[variable]/ var.t[variable])^2
    
    # nested imputation
    if(condition.m2){ 
      w.term <- 1/(m1 * (m2-1)) * ((1 - 1/m2) * var.w[variable] / var.t[variable])^2 
      lambda[variable] <- (1/m2 * var.b[variable] + (1 - 1/m2) * var.w[variable])/
        (var.u[variable] + 1/m2 * var.b[variable] + (1 - 1/m2) * var.w[variable])}
    
    # Reimputation, aopended imputation 
    else{w.term <- 0
    lambda[variable] <- (1 + 1/m1) * 1/m2 * var.b[variable] / var.t[variable]}
    
    df.complete <- n - no.parameters
    df.old[variable] <-  (b.term + w.term)^(-1)
    df.observed[variable] <- (df.complete + 1)/(df.complete + 3) * df.complete * (1 - lambda[variable])
    df[variable] <- df.old[variable] * df.observed[variable]/(df.old[variable] + df.observed[variable])
    
    if(lambda[variable] < 10^(-4)){df[variable] <- df.complete}
    }
  
  return(df)}