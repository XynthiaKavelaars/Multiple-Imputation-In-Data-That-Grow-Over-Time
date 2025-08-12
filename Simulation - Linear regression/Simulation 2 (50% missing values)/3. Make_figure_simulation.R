#### 3. Make figure simulation ####
# Make sure folder "Simulation" is set as working directory
load("Workspaces/Simulation.RData")

require(abind)
require(extrafont)#install.packages("extrafont")
font_import()
loadfonts(device = "win")
par(family = "LM Roman 10")

rho.vector <- c(0.1, 0.3, 0.5, 0.7)
#### Functions ####
results.function <- function(type, rho.vector, mp){
  results <- array(NA, dim=c(length(rho.vector), length(rho.vector),4,3))
  i <- 1
  for(rho.t1 in 1:7){
    j <- 1
    condition.rho.t1 <- (rho.t1/10) %in% rho.vector 
    if(condition.rho.t1){
    for(rho.t2 in 1:7){
      condition.rho.t2 <- (rho.t2/10) %in% rho.vector
      if(condition.rho.t2){
        filename <- paste0("out_", rho.t1, rho.t2)
        for(var in 1:4){
          for(method in 1:3){
            results[i, j, var, method] <- get(filename)[[mp]][[type]][[method]][var]
          }}
        j <- j+1}
    }
      i <- i+1}}
  dimnames(results) <- list(c("r_w 0.10", "r_w 0.30", "r_w 0.50", "r_w 0.70"),
                            c("r_b 0.10", "r_b 0.30", "r_b 0.50", "r_b 0.70"),
                            c("b0", "b_x1", "b_y1", "b_x2"),
                            c("RE-IMPUTE", "NEST", "APPEND"))
  return(results)}

  #### Extract results ####
  bias.monotone <- results.function(type="bias.lm", rho.vector=rho.vector, mp=1)
  coverage.monotone <- results.function(type="coverage.lm", rho.vector=rho.vector, mp=1)
  
  bias.nonmonotone <- results.function(type="bias.lm", rho.vector=rho.vector, mp=2)
  coverage.nonmonotone <- results.function(type="coverage.lm", rho.vector=rho.vector, mp=2)
  
  #### Print results ####
  round(bias.monotone,3)
  round(bias.nonmonotone,3)
  round(coverage.monotone,3)
  round(coverage.nonmonotone,3)
  
  #### Aggregate over similar conditions ####
  unbiased_conditions <- abind(bias.monotone, # monotone
                               bias.nonmonotone[,,,1],# nonmonotone - RE-IMPUTE
                               along=4) 
  biased_conditions <-   bias.nonmonotone[,,,2:3] # nonmonotone - NEST/APPEND
  covered_conditions <- abind(coverage.monotone, #monotone
                              coverage.nonmonotone[,,,1], # nonmonotone - RE-IMPUTE
                              along=4)
  undercovered_conditions <- coverage.nonmonotone[,,,2:3] # nonmonotone - NEST/APPEND
  
  unbiased_by1 <- apply(unbiased_conditions[,,3,], 1:2, mean)
  covered_by1 <- apply(covered_conditions[,,3,], 1:2, mean)
  biased_by1 <- apply(biased_conditions[,,3,], 1:2, mean)
  undercovered_by1 <- apply(undercovered_conditions[,,3,], 1:2, mean)
  
  #### Make figures ####
  rho.t1.vector <- rho.t2.vector <- rho.vector
  methods <- c("RE-IMPUTE", "NEST", "APPEND")
  cols<- rep("black", 4)
  lty.col <- c(1,5,3,4)
  pch.col <- c(21:24)
  
  
  setEPS()
  postscript("Plots/figure1_50.eps", family="Times")
  
  layout(matrix(c(1:9,10,10,10), nrow=4, byrow=TRUE), heights=c(2, rep(6,2),2), widths=c(5,1,5))
  par(mar=rep(0.1,4))
  plot(NULL, bty="n", xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n")
  legend("center", legend = c("Monotone (all methods)", "Nonmonotone (RE-IMPUTE)"), 
         bty = "n", cex=1.5)
  #text(x=0.5, y=0.65, label="Monotone (all methods)", cex=1.25)
  #text(x=0.5, y=0.35, label="Nonmonotone (RE-IMPUTE)", cex=1.25)
  
  plot.new()
  
  plot(NULL, bty="n", xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n")
 # text(x=0.5, y=0.5, label="Nonmonotone (NEST/APPEND)", cex=1.25)
  legend("center", legend = c("Nonmonotone (NEST/APPEND)", ""), 
         bty = "n", cex=1.5)
  
  par(mar=c(3.5,5,2,0.5))
  # Bias - unbiased conditions  
  plot(NULL, xlim=range(rho.vector), ylim=c(-0.60,0.10), 
       xlab="n", ylab="Bias",
       xaxt="n", yaxt="n", mgp=c(3.75,0.5,0), 
       main="Bias", cex.lab=1.5)
  axis(1, cex.axis=1.5, padj=-0.5, at=seq(0.1, 0.7, 0.2), cex=1.5)
  axis(2, cex.axis=1.5, las=2)
  title(xlab=expression(rho[within]), line=2, cex.lab=1.5)
  for(rho.t2 in 1:length(rho.vector)){
    points.vector <- unbiased_by1[,rho.t2]
    points(rho.t1.vector, points.vector, pch=pch.col[rho.t2], col=cols[rho.t2],
           bg=cols[rho.t2], cex=1.5)
    lines(rho.t1.vector, points.vector, col=cols[rho.t2], 
          lty=lty.col[rho.t2], lwd=2)
  }
  
  par(mar=rep(0.1,4))
  plot.new()
  
  # Bias - biased conditions  
  par(mar=c(3.5,5,2,0.5))
  plot(NULL, xlim=range(rho.vector), ylim=c(-0.60,0.10), 
       xlab="n", ylab="Bias",
       xaxt="n", yaxt="n", mgp=c(3.75,0.5,0), 
       main="Bias", cex.lab=1.5)
  axis(1, cex.axis=1.5, padj=-0.5, at=seq(0.1, 0.7, 0.2), cex=1.5)
  axis(2, cex.axis=1.5, las=2)
  title(xlab=expression(rho[within]), line=2, cex.lab=1.5)
  for(rho.t2 in 1:length(rho.vector)){
    points.vector <- biased_by1[,rho.t2]
    points(rho.t1.vector, points.vector, pch=pch.col[rho.t2], col=cols[rho.t2],
           bg=cols[rho.t2], cex=1.5)
    lines(rho.t1.vector, points.vector, col=cols[rho.t2], 
          lty=lty.col[rho.t2], lwd=2)
  }
 
  # Coverage - properly covered conditions  
  par(mar=c(3.5,5,2,0.5))
  plot(NULL, xlim=range(rho.vector), ylim=c(0,100), 
       xlab="n", ylab="Percentage",
       xaxt="n", yaxt="n", mgp=c(3.75,0.5,0), 
       main="Coverage", cex.lab=1.5)
  axis(1, cex.axis=1.5, padj=-0.5, at=seq(0.1, 0.7, 0.2), cex=1.5)
  axis(2, cex.axis=1.5, las=2)
  title(xlab=expression(rho[within]), line=2, cex.lab=1.5)
  for(rho.t2 in 1:length(rho.vector)){
    points.vector <- covered_by1[,rho.t2]
    points(rho.t1.vector, points.vector, pch=pch.col[rho.t2], col=cols[rho.t2],
           bg=cols[rho.t2], cex=1.5)
    lines(rho.t1.vector, points.vector, col=cols[rho.t2], 
          lty=lty.col[rho.t2], lwd=2)
  }
  
  par(mar=rep(0.1,4))
  plot.new()
  
  # Coverage - undercovered conditions  
  par(mar=c(3.5,5,2,0.5))
  plot(NULL, xlim=range(rho.vector), ylim=c(0,100), 
       xlab="n", ylab="Percentage",
       xaxt="n", yaxt="n", mgp=c(3.75,0.5,0), 
       main="Coverage", cex.lab=1.5)
  axis(1, cex.axis=1.5, padj=-0.5, at=seq(0.1, 0.7, 0.2), cex=1.5)
  axis(2, cex.axis=1.5, las=2)
  title(xlab=expression(rho[within]), line=2, cex.lab=1.5)
  for(rho.t2 in 1:length(rho.vector)){
    points.vector <- undercovered_by1[,rho.t2]
    points(rho.t1.vector, points.vector, pch=pch.col[rho.t2], col=cols[rho.t2],
           bg=cols[rho.t2], cex=1.5)
    lines(rho.t1.vector, points.vector, col=cols[rho.t2], 
          lty=lty.col[rho.t2], lwd=2)
  }
  
####  Legend ####

cols <- rep("black", 4)
lty.col <- c(1,5,3,4)
pch.col <- c(21:24)
rho.vector <- c(0.1, 0.3, 0.5, 0.7)
par(mar=rep(0.1,4))
plot.new()
legend("center", legend=rho.vector, col=cols, 
       lty=lty.col, title=expression(rho[between]), lwd=1,
       ncol=4, pch=pch.col, pt.bg=cols, cex=1.5)
invisible(dev.off())

