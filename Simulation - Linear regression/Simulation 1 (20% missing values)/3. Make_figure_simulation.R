#### 3. Make figure simulation ####
load("Workspaces/Simulation.RData")
require(extrafont)
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
  return(results)}

coverage.function <- function(type, rho.vector, mp){
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
              results[i, j, var, method] <- get(filename)[[mp]][[type]][var, method]
            }}
          j <- j+1}
      }
      i <- i+1}}
  return(results)}

plot.function <- function(results, type, ylim, ylab, abline=NULL, rho.vector, mp, rc, legend=TRUE){
  rho.t1.vector <- rho.t2.vector <- rho.vector
  methods <- c("reimp", "nmi", "smi")
  #cols=c("black", "gray20", "gray40", "gray60")
  cols=c("blue", "red", "green4", "magenta")
  cols<- rep("black", 4)
  lty.col <- c(1,5,3,4)
  pch.col <- c(21:24)
  mains <- c(expression(bold(b[0])), expression(bold(b[x[1]])), expression(bold(b[y[1]])),
             expression(bold(b[x[2]])))
  for(method in 2){
    filename <- paste0("Plots/", type,".png")
    png(filename = filename, family="Serif", width=600, height = 600, res = 150)
    par(mar=c(3.5,5,2,0.5))
 if(legend)   layout(c(1,2), heights=c(7,2))
    for(var in rc){
      plot(NULL, xlim=c(0, 0.7), ylim=ylim, xlab="n", ylab=ylab,
           xaxt="n", yaxt="n", mgp=c(3.75,0.5,0),            cex.main = 1.5, cex.lab=1.5)
      axis(1, cex.axis=1.5, tick=F, padj=-0.75, at=seq(0.1, 0.7, 0.2))
      axis(2, cex.axis=1.5, las=2)
      title(xlab="Cross-correlation"#expression(rho[within])
            , line=2, cex.lab=1.5)
      for(rho.t2 in 1:length(rho.vector)){
        points.vector <- results[, rho.t2, var, method]
        points(rho.t1.vector, points.vector, pch=pch.col[rho.t2], col=cols[rho.t2],
               bg=cols[rho.t2], cex=1.5)
        lines(rho.t1.vector, points.vector, col=cols[rho.t2], 
              lty=lty.col[rho.t2], lwd=2)
      }
      if (!is.null(abline)) abline(h=abline, lwd=1, lty=2, col="black")
    }
if(legend){    par(mar=c(0,0,0,0))
    plot.new()
    legend("center", legend=rho.vector, col=cols, 
           lty=lty.col, title="(Cross-)autocorrelation"#expression(rho[between])
             , lwd=2,
           ncol=3, pch=pch.col, cex=1.5)}
    dev.off()
  }
}

####  Legend ####

png(filename = "Plots/Legend.png", family="Serif", width=800, height = 150, res = 150)
#cols=c("black", "gray20", "gray40", "gray60")
cols=c("blue", "red", "green4", "magenta")
cols <- rep("black", 4)
lty.col <- c(1,5,3,4)
pch.col <- c(21:24)
rho.vector <- c(0.1, 0.3, 0.5, 0.7)
par(mar=rep(0.1,4))
plot.new()
legend("center", legend=rho.vector, col=cols, 
       lty=lty.col, title="(Cross-)autocorrelation"#expression(rho[between])
       , lwd=1,
       ncol=4, pch=pch.col, cex=1, pt.bg=cols)
dev.off()

#### * Dropout ####
# Bias
bias.monotone <- results.function(type="bias.lm", rho.vector=rho.vector, mp=1)
plot.function(bias.monotone, "Unbiased", ylim=c(-0.6,0.6), ylab = "Bias", abline=0,
              rho.vector=rho.vector, mp=1, rc=3, legend=FALSE)

# Coverage
coverage.monotone <- results.function(type="coverage.lm", rho.vector=rho.vector, mp=1)
plot.function(coverage.monotone, "Covered", ylim = c(0,100), ylab = "Coverage %", 
              abline=95, rho.vector=rho.vector, mp=1, rc=3, legend=FALSE)

#### * Non-dropout ####

bias.nonmonotone <- results.function(type="bias.lm", rho.vector=rho.vector, mp=2)
plot.function(bias.nonmonotone, "Biased", ylim = c(-0.6, 0.6), ylab = "Bias", abline=0, 
              rho.vector=rho.vector, mp=2, rc=3, legend=FALSE)

# Coverage
coverage.nonmonotone <- results.function(type="coverage.lm", rho.vector=rho.vector, mp=2)
plot.function(coverage.nonmonotone, "Undercovered", ylim = c(0,100), ylab = "Coverage %", 
              abline=95, rho.vector=rho.vector, mp=2, rc=3, legend=FALSE)

