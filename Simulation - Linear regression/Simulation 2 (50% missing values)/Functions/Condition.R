########################### Write conditions files ###########################
#mu <- mu
#rho <- rho
#nSim <- nSim
#seed <- set.seed(1545)
#n <- n
#MisPerc <- MisPerc
#m1 <- m1
#m2 <- m2
#maxit <- maxit

for(within in 1:length(rho)){
  for(between in 1:length(rho)){
    filename <- paste0("Conditions/r_0", abs(rho[within])*100, "_0", abs(rho[between])*100, ".R") # Create filename
    condition.txt <- paste0("# Condition: Rho t1: ", abs(rho[within]), ", rho t2: ", abs(rho[between])) # Name condition in condition file
    log.txt <- paste0("write(c(\" \", \"", condition.txt, "\"), \"mice_log.txt\", append=TRUE)") # Indicate condition name in log file
    mu.txt <- paste0("mu <- c(", toString(mu), ")")
    RhoWithin.txt <- paste0("w1 <- w2 <- ", rho[within])
    RhoBetween.txt <- paste0("w2 <- b1 <- b2 <- c1 <- c2 <- ", rho[between])
    nSim.txt <- paste0("nSim <- c(", toString(nSim), ")")
    n.txt <- paste0("n <- c(", toString(n), ")")
    m1.txt <- paste0("m1 <- c(", toString(m1), ")")
    m2.txt <- paste0("m2 <- c(", toString(m2), ")")
    maxit.txt <- paste0("maxit <- c(", toString(maxit), ")")
    MisPerc.txt <- paste0("MisPerc <- c(", toString(MisPerc), ")")
    sigma.txt <- paste0("sigma <- sigma.function(w1, w2, b1, b2, c1, c2)")
    resample.txt <- paste0("out_", abs(rho[within])*10, abs(rho[between])*10, "<- resample(mu, sigma, n = n, nSim = nSim, MisPerc=MisPerc, m1=m1, m2 = m2, maxit = maxit)")
    #workspace.txt <- paste0("save.image(\"Workspaces/Simulation.RData\")")
    write.vector <- c(condition.txt, " ", log.txt, " ", mu.txt, " ", RhoWithin.txt, RhoBetween.txt, " ", sigma.txt, " ", 
                      nSim.txt, n.txt, m1.txt, m2.txt, maxit.txt, MisPerc.txt, 
                      " ", resample.txt, " ")
                      #, workspace.txt)
    write(write.vector, filename)
    }
}

