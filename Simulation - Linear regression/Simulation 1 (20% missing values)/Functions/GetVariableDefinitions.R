#### variable definitions ####
mu      <- rep(0, 4)   # Define standard multivariate normal distribution 
rho     <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70) # Vary correlation structure of data
seed    <- 1545
n       <- 425 # Sample size
nSim    <- 2000 # Number of simulations 
MisPerc <- 0.2 # Missingness percentage
m1      <- 5 # Number of imputations stage 1
m2      <- 5 # Number of imputations stage 2
maxit   <- 25 # Number of iterations MICE


# Specify missingness patterns
mis.pat <- vector("list", 2)

# Monotone
mis.pat[[1]] <- matrix(c(1, 1, 1, 0,
                         1, 1, 0, 1,
                         1, 1, 0, 0,
                         1, 0, 0, 0,
                         0, 1, 0, 0),
                       nrow=5, ncol=4, byrow=TRUE)
# Non-monotone
mis.pat[[2]] <- matrix(c(1, 1, 1, 0,
                         1, 1, 0, 1,
                         1, 1, 0, 0,
                         1, 0, 1, 1,
                         1, 0, 1, 0,
                         1, 0, 0, 1,
                         1, 0, 0, 0,
                         0, 1, 1, 1,
                         0, 1, 1, 0,
                         0, 1, 0, 1,
                         0, 1, 0, 0),
                       nrow=11, ncol=4, byrow=TRUE)
pred <- matrix(c(rep(0,6),
                 rep(0,6),
                 rep(0,6),
                 rep(0,6),
                 c(0,0,1,1,0,1),
                 c(0,0,1,1,1,0)),
                 nrow=6, ncol=6, byrow=TRUE)
