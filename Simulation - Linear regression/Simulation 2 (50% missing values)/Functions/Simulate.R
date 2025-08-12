simulate.data <- function(mu, sigma, n){
  data <- as.data.frame(mvrnorm(n = n, mu = mu, Sigma = sigma, empirical = FALSE))
  colnames(data) <- c("x1", "y1", "x2", "y2")
  assign("data", data, envir = .GlobalEnv)}

