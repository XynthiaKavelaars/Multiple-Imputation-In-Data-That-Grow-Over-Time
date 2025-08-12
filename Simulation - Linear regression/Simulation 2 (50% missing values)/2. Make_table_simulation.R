#### 2. Make table simulation ####
# Make sure folder "Simulation" is set as working directory
rm(list=ls())
setwd("C:/Users/kavelaar/surfdrive/Werk/Various/Other publications/Thesis/Multivariate Behavioral Research/Working folder/Simulation 1 (50% missing values)")
#setwd("E:/Users/kavelaar/Documents/Blade/Thesis - MBR/Simulation 2 (50% missing values)")

# Load workspace
load("Workspaces/Simulation.RData")

# Load packages
require(extrafont)
require(xtable)

# Load fonts
loadfonts(device = "win")
par(family = "LM Roman 10")


rho.vector <- c(0.1, 0.3, 0.5, 0.7)
#### 2.0.1 Functions ####
results.function <- function(type, rho.vector, mp, rc){
  pop.lm.mat <- matrix(NA, nrow=length(rho.vector)*4, ncol=length(rho.vector))
  population.mu <- rep(0, 4)
  i <- 1
  for(rho.t1 in 1:7){
    j <- 1
    condition.rho.t1 <- (rho.t1/10) %in% rho.vector 
    if(condition.rho.t1){
      for(rho.t2 in 1:7){
        condition.rho.t2 <- (rho.t2/10) %in% rho.vector
        if(condition.rho.t2){
          if(type== "bias.lm"){if(rho.t2 == 7){rho.t2.i <- 6.6}
            else{rho.t2.i <- rho.t2}
            sigma <- matrix(c(10, rho.t1, rho.t2, rho.t2.i,
                              rho.t1, 10, rho.t2.i, rho.t2,
                              rho.t2, rho.t2.i, 10, rho.t2,
                              rho.t2.i, rho.t2, rho.t2, 10),
                            nrow=4, ncol=4, byrow=TRUE)
            sigma <- sigma/10
            cov.x <- sigma[c(1:3), c(1:3)]
            cov.xy <- sigma[4, c(1:3)]
            beta.hat <- solve(cov.x, cov.xy)
            intercept <- population.mu[4] - population.mu[c(1:3)] %*% beta.hat
            population.lm <- c(intercept, beta.hat)}
          else {population.lm <- rep(0, 4)}
          
          for(var in rc){
            if(type=="bias.lm"){pop.lm.mat[(i-1)*4+j, 1] <- population.lm[var]}
            filename <- paste0("out_", rho.t1, rho.t2)
          
            for(method in 1:3){
              if(type=="bias.lm"){
                pop.lm.mat[(i-1)*4+j, method+1] <- population.lm[var] + get(filename)[[mp]][[type]][[method]][var]}
              else{pop.lm.mat[(i-1)*4+j, method+1] <- get(filename)[[mp]][[type]][[method]][var]}
            }
            }
          j <- j+1}
      }
      i <- i+1}}
pop.lm.mat <- round(pop.lm.mat, 2)
rownames(pop.lm.mat) <- rep(c("$\\rho=0.1$", "$\\rho=0.3$", "$\\rho=0.5$", "$\\rho=0.7$"), length(rho.vector))
colnames(pop.lm.mat) <- c(" ", "RMI", "NMI", "AMI")
return(pop.lm.mat)}

#### 2.1 Bias, coverage, CI witdh ####
#### 2.1.1 Table settings ####

align.vector <- c("llrrrrp{0.02\\textwidth}rrrp{0.02\\textwidth}rrr")
addtorow <- list()
addtorow$pos <- list(0, 4, 8, 12)

addtorow$command <- c("\\multicolumn{1}{l}{Scenario} & \\multicolumn{4}{l}{Regression coefficients} & & \\multicolumn{3}{l}{Coverage 95\\% CI} & & \\multicolumn{3}{l}{Width 95\\% CI}\\\\ \\hline                               
                               $\\rho_{within}=0.1$ & Pop & RMI & NMI & AMI & & RMI & NMI & AMI & & RMI & NMI & AMI\\\\ \n",
                               "\\hline 
                               $\\rho_{within}=0.3$ & Pop & RMI & NMI & AMI & & RMI & NMI & AMI & & RMI & NMI & AMI\\\\ \\hline \n",
                               "\\hline
                               $\\rho_{within}=0.5$ & Pop & RMI & NMI & AMI & & RMI & NMI & AMI & & RMI & NMI & AMI\\\\ \\hline \n",
                               "\\hline
                               $\\rho_{within}=0.7$ & Pop & RMI & NMI & AMI & & RMI & NMI & AMI & & RMI & NMI & AMI\\\\ \\hline \n")
                              

row.names <- cbind(c("1.","2.","3.","4.",
                     "5.","6.","7.","8.",
                     "9.","10.","11.","12.",
                     "13.","14.","15.","16."),
                   rep(c("$\\rho_{between}=0.1$", "$\\rho_{between}=0.3$", 
                   "$\\rho_{between}=0.5$", "$\\rho_{between}=0.7$"), 4))

#### 2.1.2 Results Bias ####

abs.bias.mv.nonmonotone.by1 <- results.function(type = "bias.lm", rho.vector = rho.vector, mp = 2, rc=3)
abs.bias.mv.nonmonotone.by1 <- cbind(row.names, formatC(abs.bias.mv.nonmonotone.by1, digits = 2, format = "f"))
abs.bias.mv.monotone.by1 <- results.function(type = "bias.lm", rho.vector = rho.vector, mp = 1, rc=3)
abs.bias.mv.monotone.by1 <- cbind(row.names, formatC(abs.bias.mv.monotone.by1, digits = 2, format = "f"))




#### 2.1.3 Results Coverage ####
coverage.mv.nonmonotone.by1 <- results.function(type = "coverage.lm", rho.vector = rho.vector, mp = 2, rc=3)
coverage.mv.nonmonotone.by1 <- cbind(NA, formatC(coverage.mv.nonmonotone.by1[,2:4], digits = 1, format = "f"))
coverage.mv.monotone.by1 <- results.function(type = "coverage.lm", rho.vector = rho.vector, mp = 1, rc=3)
coverage.mv.monotone.by1 <- cbind(NA, formatC(coverage.mv.monotone.by1[,2:4], digits = 1, format = "f"))

#### 2.1.4 Results Width ####
width.mv.nonmonotone.by1 <- results.function(type = "width.lm", rho.vector = rho.vector, mp = 2, rc=3)
width.mv.nonmonotone.by1 <- cbind(NA, formatC(width.mv.nonmonotone.by1[,2:4], digits = 2, format = "f"))
width.mv.monotone.by1 <- results.function(type = "width.lm", rho.vector = rho.vector, mp = 1, rc=3)
width.mv.monotone.by1 <- cbind(NA, formatC(width.mv.monotone.by1[,2:4], digits = 2, format = "f"))

#### 2.1.5 Print tables ####
# Separate tables for monotone and nonmonotone
print.xtable(xtable(cbind(abs.bias.mv.nonmonotone.by1, coverage.mv.nonmonotone.by1, width.mv.nonmonotone.by1), 
                    caption="Population regression coefficients (Pop) and regression coefficients of the 95\\% confidence interval after reimputation (RMI), nested multiple 
		imputation (NMI), and appended multiple imputation (AMI) under a non-monotone missingness pattern 
             Results are presented for the regression coefficient of the incomplete 
             $t_1$ variable ($b_{y_1}$) under different correlation structures with 50\% missing values.", 
                    label="tab:mv_nonmonotone_by1_50misperc", 
                    align = align.vector.num), add.to.row = addtorow, 
             include.rownames = F, include.colnames = F,
             caption.placement="top", NA.string=getOption("xtable.NA.string", " "),
             sanitize.text.function = function(x){x}, booktabs=TRUE)


print.xtable(xtable(cbind(abs.bias.mv.monotone.by1, coverage.mv.monotone.by1, width.mv.monotone.by1), 
                    caption="Population regression coefficients (Pop) and regression coefficients of the 95\\% confidence interval after reimputation (RMI), nested multiple 
		imputation (NMI), and appended multiple imputation (AMI) under a monotone missingness pattern 
             Results are presented for the regression coefficient of the incomplete 
             $t_1$ variable ($b_{y_1}$) under different correlation structures with 50\% missing values.", 
                    label="tab:mv_monotone_by1_50misperc", 
                    align = align.vector.num), add.to.row = addtorow, 
             include.rownames = F, include.colnames = F,
             caption.placement="top", NA.string=getOption("xtable.NA.string", " "),
             sanitize.text.function = function(x){x}, booktabs=TRUE)


#### 2.2 Relative widths ####
#### 2.2.1 Variable means ####
# Monotone missingness
width.mu.monotone <- results.function(type="width.mu", rho.vector=rho.vector, mp=1, rc=3)
relative.width.mu.monotone <- apply(width.mu.monotone[,2:4], 1, function(x) c(x[2]/x[1], x[2]/x[3], x[3]/x[1]))
# 1: NMI/RMI; 2: NMI/AMI; 3: SMI/RMI
minimum.relative.width.mu.monotone <- apply(relative.width.mu.monotone, 1, min)
maximum.relative.width.mu.monotone <- apply(relative.width.mu.monotone, 1, max)

# Nonmonotone missingness
width.mu.nonmonotone <- results.function(type="width.mu", rho.vector=rho.vector, mp=2,rc=3)
relative.width.mu.nonmonotone <- apply(width.mu.nonmonotone[,2:4], 1, function(x) c(x[2]/x[1], x[2]/x[3], x[3]/x[1]))
# 1: NMI/RMI; 2: NMI/AMI; 3: SMI/RMI
minimum.relative.width.mu.nonmonotone <- apply(relative.width.mu.nonmonotone, 1, min)
maximum.relative.width.mu.nonmonotone <- apply(relative.width.mu.nonmonotone, 1, max)

nonmonotone.width.mu <- cbind(minimum.relative.width.mu.nonmonotone, maximum.relative.width.mu.nonmonotone)
monotone.width.mu <- cbind(minimum.relative.width.mu.monotone, maximum.relative.width.mu.monotone)
rownames(monotone.width.mu) <- rownames(nonmonotone.width.mu) <- 
  c("Nested vs. reimputation", "Nested vs. appended", "Appended vs. reimputation")
colnames(monotone.width.mu) <- colnames(nonmonotone.width.mu) <- c("Minimum", "Maximum")


#### 2.2.2 Regression coefficients ####
# Monotone missingness
width.lm.monotone <- results.function(type="width.lm", rho.vector=rho.vector, mp=1,rc=3)
relative.width.lm.monotone <- apply(width.lm.monotone[,2:4], 1, function(x) c(x[2]/x[1], x[2]/x[3], x[3]/x[1]))
# 1: NMI/RMI; 2: NMI/AMI; 3: SMI/RMI
minimum.relative.width.lm.monotone <- apply(relative.width.lm.monotone, 1, min)
maximum.relative.width.lm.monotone <- apply(relative.width.lm.monotone, 1, max)

# Nonmonotone missingness
width.lm.nonmonotone <- results.function(type="width.lm", rho.vector=rho.vector, mp=2,rc=3)
relative.width.lm.nonmonotone <- apply(width.lm.nonmonotone[,2:4], 1, function(x) c(x[2]/x[1], x[2]/x[3], x[3]/x[1]))
# 1: NMI/RMI; 2: NMI/AMI; 3: SMI/RMI
minimum.relative.width.lm.nonmonotone <- apply(relative.width.lm.nonmonotone, 1, min)
maximum.relative.width.lm.nonmonotone <- apply(relative.width.lm.nonmonotone, 1, max)

nonmonotone.width.lm <- cbind(minimum.relative.width.lm.nonmonotone, maximum.relative.width.lm.nonmonotone)
monotone.width.lm <- cbind(minimum.relative.width.lm.monotone, maximum.relative.width.lm.monotone)
rownames(monotone.width.lm) <- rownames(nonmonotone.width.lm) <- 
  c("Nested vs. reimputation", "Nested vs. appended", "Appended vs. reimputation")
colnames(monotone.width.lm) <- colnames(nonmonotone.width.lm) <- c("Minimum", "Maximum")

#### 2.2.3 Presented results relative width ####
monotone.width.lm
monotone.width.mu
nonmonotone.width.mu
