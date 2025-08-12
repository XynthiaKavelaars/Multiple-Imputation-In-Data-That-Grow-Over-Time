#### 1. Execute application ####
# Make sure folder "Application" is set as working directory
set.seed(13579)

# Load packages
require(mice)
require(MASS)

# Load functions
source("Functions/Df.R")

# Application parameters
m1 <- 10
m2.nmi <- 10
m2.smi <- 1

#### Data preparation ####
# Select variables of interest
pops.totaal <- read.csv("POPS_data.csv", header=TRUE)
outcome <- c("iq", "LHS_POPS19x", "Internal_Geel", "External_Geel", "l19_sd")
predictors.control <- c("mag", "mag2", "rok", "var148cat2", "sestot", "thsds", "par", "zwduur", "zwduur2",
                        "geslr", "ras")
predictors.weight <- c("G_SD_GA1yr", "G_SD_GA0yr")
predictors.length <- c("L_SD_GA1yr", "L_SD_GA0yr")
predictors.hc <- c("H_SD_GA1yr", "H_SD_GA0yr")
predictors.wl <- c("GL_SD_GA1yr", "GL_SD_GA0yr")
predictors.catchup <- rbind(predictors.weight, predictors.length, predictors.hc, predictors.wl)
predictors <- c(predictors.catchup, predictors.control)
predictor.names <- c("Weight", "Length", "Head circumference", "Weight-length")
selection.vars <- c("sgalen", "sgagew", "sgaho", "sgagewlen")
variables <- c("patid", outcome, predictors)
pops.data <- pops.totaal[,variables] 
pops.data <- apply(pops.data, 2, as.numeric)

# Separate selection variables: should not be imputed
selection.data <- pops.totaal[,selection.vars]
pops.data.select <- cbind(pops.data, selection.data)
# Prepare selection data for merging with imputed data
selection.data.long <- selection.data[rep(1:nrow(selection.data), times = m1),]

#### Correlations ####
upper.cor.pops <- cor.pops <- round(cor(pops.data[,c(predictors.catchup[,1],outcome)], use="complete.obs"), 3)
row.names.cor <- col.names.cor <- c("$Weight_{1yr}$", "$Length_{1yr}$", "$HC_{1yr}$", "$WL_{1yr}$",
                                    "Cognition", "HRQoL", "Internalizing behavior", 
                                    "Externalizing behavior", "$Length_{19yr}$")
upper.cor.pops <- formatC(cor.pops, digits=2, format="f")
upper.cor.pops[upper.tri(upper.cor.pops, diag=FALSE)]<-""
diag(upper.cor.pops) <- "1"


#### Missingness ####
miss.cases.adj <- miss.cases.unadj <- miss.values.adj <- miss.values.unadj <- 
  perc.values.congenial.list.adj <- perc.values.congenial.list.unadj <- vector("list", length(outcome))

for(out in 1:length(outcome)){
  #### Weight data ####
  sel.data.weight <- pops.data.select[!is.na(pops.data.select[,"sgagew"]) & pops.data.select[,"sgagew"] == 1, 
                               c("patid", selection.vars, predictors, outcome[out])]
  mp.weight <- as.data.frame(md.pattern(sel.data.weight))
  mp.weight[order(c(outcome[out][1:length(outcome[out])])),]
  colnames(mp.weight)[ncol(mp.weight)] <- "n.miss" 
  no.cases <- rownames(mp.weight)
  rownames(mp.weight) <- c(1:nrow(mp.weight))
  mp.weight <- cbind(mp.weight, no.cases)
  mp.weight <- apply(mp.weight, 2, as.numeric)
  
  n.miss.weight <- sum(is.na(sel.data.weight))
  n.combs.weight <- nrow(mp.weight)-1
  mp.weight[mp.weight[,outcome[out]]==0,]
  complete.weight.adj <- mp.weight[which(apply(mp.weight, 1, function(x) all(x[c(predictors.weight,predictors.control,outcome[out])] == 1))),]
  congenial.weight.adj <- mp.weight[which(mp.weight[,outcome[out]]==0),]
  uncongenial.weight.adj <- mp.weight[which(mp.weight[,outcome[out]]==1 &
                                              apply(mp.weight, 1, function(x) any(x[c(predictors.weight,predictors.control)] == 0))),]
  
  complete.weight.unadj <- mp.weight[which(apply(mp.weight, 1, function(x) all(x[c(predictors.weight,outcome[out])] == 1))),]
  congenial.weight.unadj <- mp.weight[which(mp.weight[,outcome[out]]==0),]
  uncongenial.weight.unadj <- mp.weight[which(mp.weight[,outcome[out]]==1 &
                                                apply(mp.weight, 1, function(x) any(x[c(predictors.weight)] == 0))),]
  
  n.monotone.combs.weight.adj <- nrow(congenial.weight.adj)
  n.cases.monotone.weight.adj <- sum(congenial.weight.adj[,"no.cases"])
  n.values.monotone.weight.adj <- sum(congenial.weight.adj[,"n.miss"]*congenial.weight.adj[,"no.cases"])
  
  n.monotone.combs.weight.unadj <- nrow(congenial.weight.unadj)
  n.cases.monotone.weight.unadj <- sum(congenial.weight.unadj[,"no.cases"])
  n.values.monotone.weight.unadj <- sum(congenial.weight.unadj[,"n.miss"]*congenial.weight.unadj[,"no.cases"])
  
  n.nonmonotone.combs.weight.adj <- nrow(uncongenial.weight.adj)
  n.cases.nonmonotone.weight.adj <- sum(uncongenial.weight.adj[,"no.cases"])
  n.values.nonmonotone.weight.adj <- sum(uncongenial.weight.adj[,"n.miss"]*uncongenial.weight.adj[,"no.cases"])
  
  n.nonmonotone.combs.weight.unadj <- nrow(uncongenial.weight.unadj)
  n.cases.nonmonotone.weight.unadj <- sum(uncongenial.weight.unadj[,"no.cases"])
  n.values.nonmonotone.weight.unadj <- sum(uncongenial.weight.unadj[,"n.miss"]*uncongenial.weight.unadj[,"no.cases"])
  
  n.complete.weight.adj <- nrow(complete.weight.adj)
  n.cases.complete.weight.adj <- sum(complete.weight.adj[,"no.cases"])
  
  n.complete.weight.unadj <- nrow(complete.weight.unadj)
  n.cases.complete.weight.unadj <- sum(complete.weight.unadj[,"no.cases"])
  

  #### Length data ####
  sel.data.length <- pops.data.select[!is.na(pops.data.select[,"sgalen"]) & pops.data.select[,"sgalen"] == 1,
                               c("patid", selection.vars, predictors, outcome[out])]
  mp.length <- as.data.frame(md.pattern(sel.data.length))
  mp.length[order(c(outcome[out][1:length(outcome[out])])),]
  colnames(mp.length)[ncol(mp.length)] <- "n.miss" 
  no.cases <- rownames(mp.length)
  rownames(mp.length) <- c(1:nrow(mp.length))
  mp.length <- cbind(mp.length, no.cases)
  mp.length <- apply(mp.length, 2, as.numeric)
  
  n.miss.length <- sum(is.na(sel.data.length))
  n.combs.length <- nrow(mp.length)-1
  
  complete.length.adj <- mp.length[which(apply(mp.length, 1, function(x) all(x[c(predictors.length,predictors.control,outcome[out])] == 1))),]
  congenial.length.adj <- mp.length[which(mp.length[,outcome[out]]==0),]
  uncongenial.length.adj <- mp.length[which(mp.length[,outcome[out]]==1 &
                                              apply(mp.length, 1, function(x) any(x[c(predictors.length,predictors.control)] == 0))),]
  
  complete.length.unadj <- mp.length[which(apply(mp.length, 1, function(x) all(x[c(predictors.length,outcome[out])] == 1))),]
  congenial.length.unadj <- mp.length[which(mp.length[,outcome[out]]==0),]
  uncongenial.length.unadj <- mp.length[which(mp.length[,outcome[out]]==1 &
                                                apply(mp.length, 1, function(x) any(x[c(predictors.length)] == 0))),]
  
  n.monotone.combs.length.adj <- nrow(congenial.length.adj)
  n.cases.monotone.length.adj <- sum(congenial.length.adj[,"no.cases"])
  n.values.monotone.length.adj <- sum(congenial.length.adj[,"n.miss"]*congenial.length.adj[,"no.cases"])
  
  n.monotone.combs.length.unadj <- nrow(congenial.length.unadj)
  n.cases.monotone.length.unadj <- sum(congenial.length.unadj[,"no.cases"])
  n.values.monotone.length.unadj <- sum(congenial.length.unadj[,"n.miss"]*congenial.length.unadj[,"no.cases"])
  
  n.nonmonotone.combs.length.adj <- nrow(uncongenial.length.adj)
  n.cases.nonmonotone.length.adj <- sum(uncongenial.length.adj[,"no.cases"])
  n.values.nonmonotone.length.adj <- sum(uncongenial.length.adj[,"n.miss"]*uncongenial.length.adj[,"no.cases"])
  
  n.nonmonotone.combs.length.unadj <- nrow(uncongenial.length.unadj)
  n.cases.nonmonotone.length.unadj <- sum(uncongenial.length.unadj[,"no.cases"])
  n.values.nonmonotone.length.unadj <- sum(uncongenial.length.unadj[,"n.miss"]*uncongenial.length.unadj[,"no.cases"])
  
  n.complete.length.adj <- nrow(complete.length.adj)
  n.cases.complete.length.adj <- sum(complete.length.adj[,"no.cases"])
  
  n.complete.length.unadj <- nrow(complete.length.unadj)
  n.cases.complete.length.unadj <- sum(complete.length.unadj[,"no.cases"])
  #### Head circumference data ####
  sel.data.hc <- pops.data.select[!is.na(pops.data.select[,"sgaho"]) & pops.data.select[,"sgaho"] == 1,
                           c("patid", selection.vars, predictors, outcome[out])]
  mp.hc <- as.data.frame(md.pattern(sel.data.hc))
  mp.hc[order(c(outcome[out][1:length(outcome[out])])),]
  colnames(mp.hc)[ncol(mp.hc)] <- "n.miss" 
  no.cases <- rownames(mp.hc)
  rownames(mp.hc) <- c(1:nrow(mp.hc))
  mp.hc <- cbind(mp.hc, no.cases)
  mp.hc <- apply(mp.hc, 2, as.numeric)
  
  n.miss.hc <- sum(is.na(sel.data.hc))
  n.combs.hc <- nrow(mp.hc)-1
  
  complete.hc.adj <- mp.hc[which(apply(mp.hc, 1, function(x) all(x[c(predictors.hc,predictors.control,outcome[out])] == 1))),]
  congenial.hc.adj <- mp.hc[which(mp.hc[,outcome[out]]==0),]
  uncongenial.hc.adj <- mp.hc[which(mp.hc[,outcome[out]]==1 &
                                      apply(mp.hc, 1, function(x) any(x[c(predictors.hc,predictors.control)] == 0))),]
  
  complete.hc.unadj <- mp.hc[which(apply(mp.hc, 1, function(x) all(x[c(predictors.hc,outcome[out])] == 1))),]
  congenial.hc.unadj <- mp.hc[which(mp.hc[,outcome[out]]==0),]
  uncongenial.hc.unadj <- mp.hc[which(mp.hc[,outcome[out]]==1 &
                                        apply(mp.hc, 1, function(x) any(x[c(predictors.hc)] == 0))),]
  
  n.monotone.combs.hc.adj <- nrow(congenial.hc.adj)
  n.cases.monotone.hc.adj <- sum(congenial.hc.adj[,"no.cases"])
  n.values.monotone.hc.adj <- sum(congenial.hc.adj[,"n.miss"]*congenial.hc.adj[,"no.cases"])
  
  n.monotone.combs.hc.unadj <- nrow(congenial.hc.unadj)
  n.cases.monotone.hc.unadj <- sum(congenial.hc.unadj[,"no.cases"])
  n.values.monotone.hc.unadj <- sum(congenial.hc.unadj[,"n.miss"]*congenial.hc.unadj[,"no.cases"])
  
  n.nonmonotone.combs.hc.adj <- nrow(uncongenial.hc.adj)
  n.cases.nonmonotone.hc.adj <- sum(uncongenial.hc.adj[,"no.cases"])
  n.values.nonmonotone.hc.adj <- sum(uncongenial.hc.adj[,"n.miss"]*uncongenial.hc.adj[,"no.cases"])
  
  n.nonmonotone.combs.hc.unadj <- nrow(uncongenial.hc.unadj)
  n.cases.nonmonotone.hc.unadj <- sum(uncongenial.hc.unadj[,"no.cases"])
  n.values.nonmonotone.hc.unadj <- sum(uncongenial.hc.unadj[,"n.miss"]*uncongenial.hc.unadj[,"no.cases"])
  
  n.complete.hc.adj <- nrow(complete.hc.adj)
  n.cases.complete.hc.adj <- sum(complete.hc.adj[,"no.cases"])
  
  n.complete.hc.unadj <- nrow(complete.hc.unadj)
  n.cases.complete.hc.unadj <- sum(complete.hc.unadj[,"no.cases"])
  #### weight-length data ####
  sel.data.wl <- pops.data.select[!is.na(pops.data.select[,"sgagewlen"]) & pops.data.select[,"sgagewlen"] == 1,
                           c("patid", selection.vars, predictors, outcome[out])]
  mp.wl <- as.data.frame(md.pattern(sel.data.wl))
  mp.wl[order(c(outcome[out][1:length(outcome[out])])),]
  colnames(mp.wl)[ncol(mp.wl)] <- "n.miss" 
  no.cases <- rownames(mp.wl)
  rownames(mp.wl) <- c(1:nrow(mp.wl))
  mp.wl <- cbind(mp.wl, no.cases)
  mp.wl <- apply(mp.wl, 2, as.numeric)
  
  n.miss.wl <- sum(is.na(sel.data.wl))
  n.combs.wl <- nrow(mp.wl)-1
  
  complete.wl.adj <- mp.wl[which(apply(mp.wl, 1, function(x) all(x[c(predictors.wl,predictors.control,outcome[out])] == 1))),]
  congenial.wl.adj <- mp.wl[which(mp.wl[,outcome[out]]==0),]
  uncongenial.wl.adj <- mp.wl[which(mp.wl[,outcome[out]]==1 &
                                      apply(mp.wl, 1, function(x) any(x[c(predictors.wl,predictors.control)] == 0))),]
  
  complete.wl.unadj <- mp.wl[which(apply(mp.wl, 1, function(x) all(x[c(predictors.wl,outcome[out])] == 1))),]
  congenial.wl.unadj <- mp.wl[which(mp.wl[,outcome[out]]==0),]
  uncongenial.wl.unadj <- mp.wl[which(mp.wl[,outcome[out]]==1 &
                                        apply(mp.wl, 1, function(x) any(x[c(predictors.wl)] == 0))),]
  
  n.monotone.combs.wl.adj <- nrow(congenial.wl.adj)
  n.cases.monotone.wl.adj <- sum(congenial.wl.adj[,"no.cases"])
  n.values.monotone.wl.adj <- sum(congenial.wl.adj[,"n.miss"]*congenial.wl.adj[,"no.cases"])
  
  n.monotone.combs.wl.unadj <- nrow(congenial.wl.unadj)
  n.cases.monotone.wl.unadj <- sum(congenial.wl.unadj[,"no.cases"])
  n.values.monotone.wl.unadj <- sum(congenial.wl.unadj[,"n.miss"]*congenial.wl.unadj[,"no.cases"])
  
  n.nonmonotone.combs.wl.adj <- nrow(uncongenial.wl.adj)
  n.cases.nonmonotone.wl.adj <- sum(uncongenial.wl.adj[,"no.cases"])
  n.values.nonmonotone.wl.adj <- sum(uncongenial.wl.adj[,"n.miss"]*uncongenial.wl.adj[,"no.cases"])
  
  n.nonmonotone.combs.wl.unadj <- nrow(uncongenial.wl.unadj)
  n.cases.nonmonotone.wl.unadj <- sum(uncongenial.wl.unadj[,"no.cases"])
  n.values.nonmonotone.wl.unadj <- sum(uncongenial.wl.unadj[,"n.miss"]*uncongenial.wl.unadj[,"no.cases"])
  
  n.complete.wl.adj <- nrow(complete.wl.adj)
  n.cases.complete.wl.adj <- sum(complete.wl.adj[,"no.cases"])
  
  n.complete.wl.unadj <- nrow(complete.wl.unadj)
  n.cases.complete.wl.unadj <- sum(complete.wl.unadj[,"no.cases"])
  #### Make table ####
  n.complete.adj <- c(n.cases.complete.weight.adj, n.cases.complete.length.adj, n.cases.complete.hc.adj, n.cases.complete.wl.adj)
  n.total <- c(nrow(sel.data.weight), nrow(sel.data.length), nrow(sel.data.hc), nrow(sel.data.wl))
  n.complete.unadj <- c(n.cases.complete.weight.unadj, n.cases.complete.length.unadj, n.cases.complete.hc.unadj, n.cases.complete.wl.unadj)
  n.values.congenial.adj <- c(n.values.monotone.weight.adj, n.values.monotone.length.adj, n.values.monotone.hc.adj, n.values.monotone.wl.adj)
  n.values.uncongenial.adj <- c(n.values.nonmonotone.weight.adj, n.values.nonmonotone.length.adj, n.values.nonmonotone.hc.adj, n.values.nonmonotone.wl.adj)
  n.values.congenial.unadj <- c(n.values.monotone.weight.unadj, n.values.monotone.length.unadj, n.values.monotone.hc.unadj, n.values.monotone.wl.unadj)
  n.values.uncongenial.unadj <- c(n.values.nonmonotone.weight.unadj, n.values.nonmonotone.length.unadj, n.values.nonmonotone.hc.unadj, n.values.nonmonotone.wl.unadj)
  n.values.adj <- n.values.congenial.adj + n.values.uncongenial.adj
  n.values.unadj <- n.values.congenial.unadj + n.values.uncongenial.unadj

  perc.values.congenial.list.adj[[out]] <- n.values.congenial.adj/n.values.adj*100
  perc.values.congenial.list.unadj[[out]] <- n.values.congenial.unadj/n.values.unadj*100
}
n.total <- c(nrow(sel.data.weight), nrow(sel.data.length), nrow(sel.data.hc), nrow(sel.data.wl))
n.values.total.adj <- n.total*(length(predictors.control)+3)
n.values.total.unadj <- n.total*3
perc.total.adj <- n.values.adj/n.values.total.adj*100
perc.total.unadj <- n.values.unadj/n.values.total.unadj*100




miss.df.adj <- cbind(c("Weight", "Length", "HC", "WL"),
                         formatC(n.total, digits = 0, format = "f"),
                         formatC(n.complete.adj, digits = 0, format = "f"),
                         formatC(perc.total.adj, digits = 1, format = "f"),
                         formatC(perc.values.congenial.list.adj[[1]], digits = 1, format = "f"),
                         formatC(perc.values.congenial.list.adj[[2]], digits = 1, format = "f"),
                         formatC(perc.values.congenial.list.adj[[3]], digits = 1, format = "f"),
                         formatC(perc.values.congenial.list.adj[[4]], digits = 1, format = "f"),
                         formatC(perc.values.congenial.list.adj[[5]], digits = 1, format = "f"))
colnames(miss.df.adj) <- c(" ", "Total n", "Complete n", "% missing", 
                               outcome[1], outcome[2], outcome[3], outcome[4], outcome[5])

miss.df.unadj <- cbind(c("Weight", "Length", "HC", "WL"),
                           formatC(n.total, digits = 0, format = "f"),
                           formatC(n.complete.unadj, digits = 0, format = "f"),
                           formatC(perc.total.unadj, digits = 1, format = "f"),
                           formatC(perc.values.congenial.list.unadj[[1]], digits = 1, format = "f"),
                           formatC(perc.values.congenial.list.unadj[[2]], digits = 1, format = "f"),
                           formatC(perc.values.congenial.list.unadj[[3]], digits = 1, format = "f"),
                           formatC(perc.values.congenial.list.unadj[[4]], digits = 1, format = "f"),
                           formatC(perc.values.congenial.list.unadj[[5]], digits = 1, format = "f"))
colnames(miss.df.unadj) <- c(" ", "Total n", "Complete n", "% missing", 
                                 outcome[1], outcome[2], outcome[3], outcome[4], outcome[5])


#### Imputation ####

# t1.data
pops.t1 <- pops.data.select[,c("patid", predictors)]
cor.t1 <- cor(pops.t1, use="complete.obs")
min.cor.t1 <- min(abs(cor.t1))
mean.cor.t1 <- mean(cor.t1)

# t2.data
pops.t2 <- pops.data.select[,outcome]
pops.data.select <- cbind(pops.t1, pops.t2)

# Predictor matrix reimputation
pred <- matrix(rep(0, ncol(pops.data.select)^2), nrow=ncol(pops.data.select), ncol=ncol(pops.data.select), byrow=TRUE)
colnames(pred) <- rownames(pred) <- colnames(pops.data.select)

# Predictor matrix t1.imputation
pred.t1 <- matrix(rep(0, ncol(pops.t1)^2), nrow=ncol(pops.t1), ncol=ncol(pops.t1), byrow=TRUE)
colnames(pred.t1) <- rownames(pred.t1) <- colnames(pops.t1)

# Predictor matrix t2.imputation
pred.t1[predictors, predictors] <- 
  pred[c(predictors,outcome), c(predictors, outcome)] <- 1
pred.t1[c(predictors.weight, predictors.length), predictors.wl] <- 
  pred[c(predictors.weight, predictors.length), predictors.wl] <- 0
pred.t1[c(predictors.weight, predictors.length), predictors.wl[2]] <- 
  pred[c(predictors.weight, predictors.length), predictors.wl[2]] <- 0
diag(pred.t1) <- diag(pred) <- 0

pred.t2 <- rbind(rep(0, ncol(pred)+2), rep(0, ncol(pred)+2), 
                 cbind(rep(0, nrow(pred)), rep(0, nrow(pred)), pred))
colnames(pred.t2) <- rownames(pred.t2) <- c(".imp", ".id", colnames(pred))

#### Formulas####
formulas.unadjusted <- formulas.adjusted <- rep(list(list()), length(outcome))
for(out in 1:length(outcome)){
for(predictor in 1:nrow(predictors.catchup)){
  formulas.unadjusted[[out]][[predictor]] <- as.formula(paste(outcome[[out]], paste(predictors.catchup[predictor,], collapse=" + "), sep=" ~ "))
  formulas.adjusted[[out]][[predictor]] <- as.formula(paste(outcome[[out]], paste(c(predictors.catchup[predictor,], predictors.control), collapse=" + "), sep=" ~ "))}}


#### Storage objects ####
fit.control.nmi <- q.hat.nest.control.nmi <- q.hat.control.nmi <- q.nest.control.nmi <- e.w.control.nmi <-  
  w.nest.control.nmi <- u.control.nmi <-  fit.control.smi <- q.hat.control.smi <- u.control.smi <- 
  rep(list(rep(list(list()), nrow(predictors.catchup))), length(outcome)) 
q.bar.control.reimp <- df.control.reimp <- t.value.control.reimp <- p.value.control.reimp <- ci.95.control.reimp <- width.control.reimp <-
  q.bar.control.nmi <- var.u.control.nmi <- e.b.control.nmi <- var.b.control.nmi <- var.w.control.nmi <- var.t.control.nmi <- 
  df.control.nmi <- q.t.control.nmi <- t.value.control.nmi <- p.value.control.nmi <- ci.95.control.nmi <- width.control.nmi <-
  q.bar.control.smi <- var.u.control.smi <- e.b.control.smi <- var.b.control.smi <- var.w.control.smi <- var.t.control.smi <- 
  df.control.smi <- q.t.control.smi <- t.value.control.smi <- p.value.control.smi <- ci.95.control.smi <- width.control.smi <- 
  rep(list(list()), length(outcome)) 
fit.nmi <- q.hat.nest.nmi <- q.hat.nmi <- q.nest.nmi <- e.w.nmi <-  
  w.nest.nmi <- u.nmi <-  fit.smi <- q.hat.smi <- u.smi <- 
  rep(list(rep(list(list()), nrow(predictors.catchup))), length(outcome)) 
q.bar.reimp <- df.reimp <- t.value.reimp <- p.value.reimp <- ci.95.reimp <- width.reimp <-
  q.bar.nmi <- var.u.nmi <- e.b.nmi <- var.b.nmi <- var.w.nmi <- var.t.nmi <- 
  df.nmi <- q.t.nmi <- t.value.nmi <- p.value.nmi <- ci.95.nmi <- width.nmi <-
  q.bar.smi <- var.u.smi <- e.b.smi <- var.b.smi <- var.w.smi <- var.t.smi <- 
  df.smi <- q.t.smi <- t.value.smi <- p.value.smi <- ci.95.smi <- width.smi <-rep(list(list()), length(outcome)) 

#### Imputation reimputation/t1 ####
t2.reimp <- mice(pops.data, predictorMatrix=pred, maxit=20, m=m1, print=FALSE)
t2.reimp.data <- cbind(complete(t2.reimp, action="long", include=FALSE), selection.data.long)
t2.reimp.data[,selection.vars] <- apply(t2.reimp.data[,selection.vars], 2, as.numeric)
par(mfrow=c(1,1))

# Check convergence
plot(t2.reimp, layout=c(6,6))

# Imputation - NMI
t1.imp <- mice(pops.t1, predictorMatrix=pred.t1, maxit=20, m=m1, print=FALSE)
t1.imp.data <- complete(t1.imp, action="long", include=TRUE)
t2.data.long <- pops.t2[rep(1:nrow(pops.t2), times = m1+1),]
t1.t2.data <- cbind(t1.imp.data, t2.data.long)
t1.t2.data[,".imp"] <- as.numeric(as.character(t1.t2.data[,".imp"]))
colnames(t1.t2.data)[(ncol(t1.t2.data)-length(outcome)+1):ncol(t1.t2.data)] <- outcome

#### Imputation t2 ####
for(imp.m1 in 1:m1){
  pops.data.nest <- t1.t2.data[t1.t2.data$.imp == imp.m1,]
  t2.imp.nmi <- mice(pops.data.nest, predictorMatrix=pred.t2, maxit=20, m=m2.nmi, print=FALSE) 
  t2.imp.smi <- mice(pops.data.nest, predictorMatrix=pred.t2, maxit=20, m=m2.smi, print=FALSE)
  t2.nmi.data.nest <- complete(t2.imp.nmi, action="long")
  t2.smi.data.nest <- complete(t2.imp.smi, action="long")
  #### Analysis per nest: NMI ####
  for(imp.m2 in 1:m2.nmi){
    imp.data.nmi <- cbind(t2.nmi.data.nest[t2.nmi.data.nest$.imp == imp.m2,], selection.data)
    for(out in 1:length(outcome)){
      for (predictor in 1:nrow(predictors.catchup)){
      imp.data.nest.nmi <- imp.data.nmi[!is.na(imp.data.nmi[,selection.vars[predictor]]) & imp.data.nmi[,selection.vars[predictor]] == 1,]
      fit.control.nmi <- lm(formulas.adjusted[[out]][[predictor]], data = imp.data.nest.nmi)
      fit.nmi <- lm(formulas.unadjusted[[out]][[predictor]], data = imp.data.nest.nmi)
      q.hat.nest.control.nmi[[out]][[predictor]][[imp.m2]] <- fit.control.nmi$coefficients
      q.hat.nest.nmi[[out]][[predictor]][[imp.m2]] <- fit.nmi$coefficients
      u.control.nmi[[out]][[predictor]][[(imp.m1-1)*m2.nmi+imp.m2]] <- vcov(fit.control.nmi)
      u.nmi[[out]][[predictor]][[(imp.m1-1)*m2.nmi+imp.m2]] <- vcov(fit.nmi)
    }
    }}
  
  #### Analysis per nest: SMI ####
  imp.data.smi <- cbind(t2.smi.data.nest, selection.data)
  for(out in 1:length(outcome)){
  for (predictor in 1:nrow(predictors.catchup)){
    imp.data.nest.smi <- imp.data.smi[!is.na(imp.data.smi[,selection.vars[predictor]]) & imp.data.smi[,selection.vars[predictor]] == 1,]
    
    fit.control.smi <- lm(formulas.adjusted[[out]][[predictor]], data = imp.data.nest.smi)
    fit.smi <- lm(formulas.unadjusted[[out]][[predictor]], data = imp.data.nest.smi)
    q.hat.control.smi[[out]][[predictor]][[imp.m1]] <- fit.control.smi$coefficients
    q.hat.smi[[out]][[predictor]][[imp.m1]] <- fit.smi$coefficients
    u.control.smi[[out]][[predictor]][[imp.m1]] <- vcov(fit.control.smi)
    u.smi[[out]][[predictor]][[imp.m1]] <- vcov(fit.smi)
    
    ### Pool per nest: NMI ####
    q.hat.control.nmi[[out]][[predictor]][[imp.m1]] <- q.hat.nest.control.nmi[[out]][[predictor]]
    q.hat.nmi[[out]][[predictor]][[imp.m1]] <- q.hat.nest.nmi[[out]][[predictor]]
    q.nest.control.nmi[[out]][[predictor]][[imp.m1]] <- apply(simplify2array(q.hat.nest.control.nmi[[out]][[predictor]]), 1, mean)
    q.nest.nmi[[out]][[predictor]][[imp.m1]] <- apply(simplify2array(q.hat.nest.nmi[[out]][[predictor]]), 1, mean)
    e.w.control.nmi[[out]][[predictor]] <- t(simplify2array(lapply(q.hat.nest.control.nmi[[out]][[predictor]], function(x) x - q.nest.control.nmi[[out]][[predictor]][[imp.m1]])))
    e.w.nmi[[out]][[predictor]] <- t(simplify2array(lapply(q.hat.nest.nmi[[out]][[predictor]], function(x) x - q.nest.nmi[[out]][[predictor]][[imp.m1]])))
    w.nest.control.nmi[[out]][[predictor]][[imp.m1]] <- t(e.w.control.nmi[[out]][[predictor]]) %*% e.w.control.nmi[[out]][[predictor]]
    w.nest.nmi[[out]][[predictor]][[imp.m1]] <- t(e.w.nmi[[out]][[predictor]]) %*% e.w.nmi[[out]][[predictor]]
  }}
  }


#### Reimputation & pooling NMI/SMI####
fit.control.reimp <- vector("list", m1)
fit.reimp <- vector("list", m1)
for(out in 1:length(outcome)){
for(predictor in 1:nrow(predictors.catchup)){
  tryCatch({
  for(imp.m1 in 1:m1){
    imp.data.imp <- t2.reimp.data[t2.reimp.data$.imp == imp.m1,]
    imp.data.reimp <- imp.data.imp[imp.data.imp[,selection.vars[predictor]] == 1 &
                                   !is.na(imp.data.imp[,selection.vars[predictor]]) ,]
    fit.control.reimp[[imp.m1]] <- lm(formulas.adjusted[[out]][[predictor]], data = imp.data.reimp)
    fit.reimp[[imp.m1]] <- lm(formulas.unadjusted[[out]][[predictor]], data = imp.data.reimp)}
  fit.control.reimp.mira <- as.mira(fit.control.reimp)
  fit.reimp.mira <- as.mira(fit.reimp)
  pool.control.reimp <- pool(fit.control.reimp.mira)
  pool.reimp <- pool(fit.reimp.mira)
  q.bar.control.reimp[[out]][[predictor]] <- summary(pool.control.reimp)[,"est"]
  q.bar.reimp[[out]][[predictor]] <- summary(pool.reimp)[,"est"]
  ci.95.control.reimp[[out]][[predictor]] <- cbind(summary(pool.control.reimp)[,"lo 95"], 
                                             summary(pool.control.reimp)[,"hi 95"])
  ci.95.reimp[[out]][[predictor]] <- cbind(summary(pool.reimp)[,"lo 95"], 
                               summary(pool.reimp)[,"hi 95"])
  t.value.control.reimp[[out]][[predictor]] <- summary(pool.control.reimp)[,"t"]
  t.value.reimp[[out]][[predictor]] <- summary(pool.reimp)[,"t"]
  p.value.control.reimp[[out]][[predictor]] <- summary(pool.control.reimp)[,"Pr(>|t|)"]
  p.value.reimp[[out]][[predictor]] <- summary(pool.reimp)[,"Pr(>|t|)"]
  df.control.reimp[[out]][[predictor]] <- summary(pool.control.reimp)[,"df"]
  df.reimp[[out]][[predictor]] <- summary(pool.reimp)[,"df"]
  width.control.reimp[[out]][[predictor]] <- ci.95.control.reimp[[out]][[predictor]][,2] - ci.95.control.reimp[[out]][[predictor]][,1]
  width.reimp[[out]][[predictor]] <- ci.95.reimp[[out]][[predictor]][,2] - ci.95.reimp[[out]][[predictor]][,1]
  
  q.bar.control.nmi[[out]][[predictor]] <- apply(simplify2array(q.nest.control.nmi[[out]][[predictor]]), 1, mean)
  q.bar.nmi[[out]][[predictor]] <- apply(simplify2array(q.nest.nmi[[out]][[predictor]]), 1, mean)
  var.u.control.nmi[[out]][[predictor]] <- diag(apply(simplify2array(u.control.nmi[[out]][[predictor]]), c(1:2), mean))
  var.u.nmi[[out]][[predictor]] <- diag(apply(simplify2array(u.nmi[[out]][[predictor]]), c(1:2), mean))
  e.b.control.nmi[[out]][[predictor]] <- t(simplify2array(lapply(q.nest.control.nmi[[out]][[predictor]], function(x) x - q.bar.control.nmi[[out]][[predictor]])))
  e.b.nmi[[out]][[predictor]] <- t(simplify2array(lapply(q.nest.nmi[[out]][[predictor]], function(x) x - q.bar.nmi[[out]][[predictor]])))
  var.b.control.nmi[[out]][[predictor]] <- diag(t(e.b.control.nmi[[out]][[predictor]])%*%e.b.control.nmi[[out]][[predictor]] * m2.nmi/(m1-1))
  var.b.nmi[[out]][[predictor]] <- diag(t(e.b.nmi[[out]][[predictor]])%*%e.b.nmi[[out]][[predictor]] * m2.nmi/(m1-1))
  var.w.control.nmi[[out]][[predictor]] <- 1/(m1*(m2.nmi-1)) * diag(apply(simplify2array(w.nest.control.nmi[[out]][[predictor]]), c(1:2), mean))
  var.w.nmi[[out]][[predictor]] <- 1/(m1*(m2.nmi-1)) * diag(apply(simplify2array(w.nest.nmi[[out]][[predictor]]), c(1:2), mean))
  var.t.control.nmi[[out]][[predictor]] <- var.u.control.nmi[[out]][[predictor]] + 1/m2.nmi * (1 + 1/m1) * var.b.control.nmi[[out]][[predictor]] + (1 - 1/m2.nmi) * var.w.control.nmi[[out]][[predictor]]
  var.t.nmi[[out]][[predictor]] <- var.u.nmi[[out]][[predictor]] + 1/m2.nmi * (1 + 1/m1) * var.b.nmi[[out]][[predictor]] + (1 - 1/m2.nmi) * var.w.nmi[[out]][[predictor]]
  df.control.nmi[[out]][[predictor]] <- df.function(var.u.control.nmi[[out]][[predictor]], var.b.control.nmi[[out]][[predictor]], var.w.control.nmi[[out]][[predictor]], var.t.control.nmi[[out]][[predictor]], m1=m1, m2=m2.nmi, 
                                              nrow(imp.data.nest.nmi), length(q.bar.control.nmi[[out]][[predictor]]), length(q.bar.control.nmi[[out]][[predictor]]))
  df.nmi[[out]][[predictor]] <- df.function(var.u.nmi[[out]][[predictor]], var.b.nmi[[out]][[predictor]], var.w.nmi[[out]][[predictor]], var.t.nmi[[out]][[predictor]], m1=m1, m2=m2.nmi, 
                                        nrow(imp.data.nest.nmi), length(q.bar.nmi[[out]][[predictor]]), length(q.bar.nmi[[out]][[predictor]]))
  q.t.control.nmi[[out]][[predictor]] <- qt(0.975, df.control.nmi[[out]][[predictor]])
  q.t.nmi[[out]][[predictor]] <- qt(0.975, df.nmi[[out]][[predictor]])
  t.value.control.nmi[[out]][[predictor]] <- (q.bar.control.nmi[[out]][[predictor]] - 0)/sqrt(var.t.control.nmi[[out]][[predictor]])
  t.value.nmi[[out]][[predictor]] <- (q.bar.nmi[[out]][[predictor]] - 0)/sqrt(var.t.nmi[[out]][[predictor]])
  p.value.control.nmi[[out]][[predictor]] <- round(pt(t.value.control.nmi[[out]][[predictor]], df.control.nmi[[out]][[predictor]], lower.tail=FALSE), 3)
  p.value.nmi[[out]][[predictor]] <- round(pt(t.value.nmi[[out]][[predictor]], df.nmi[[out]][[predictor]], lower.tail=FALSE), 3)
  ci.95.control.nmi[[out]][[predictor]] <- round(cbind(q.bar.control.nmi[[out]][[predictor]] - q.t.control.nmi[[out]][[predictor]] * sqrt(var.t.control.nmi[[out]][[predictor]]), 
                                                 q.bar.control.nmi[[out]][[predictor]] + q.t.control.nmi[[out]][[predictor]] * sqrt(var.t.control.nmi[[out]][[predictor]])), 3)
  ci.95.nmi[[out]][[predictor]] <- round(cbind(q.bar.nmi[[out]][[predictor]] - q.t.nmi[[out]][[predictor]] * sqrt(var.t.nmi[[out]][[predictor]]), 
                                           q.bar.nmi[[out]][[predictor]] + q.t.nmi[[out]][[predictor]] * sqrt(var.t.nmi[[out]][[predictor]])), 3)
  width.control.nmi[[out]][[predictor]] <- ci.95.control.nmi[[out]][[predictor]][,2] - ci.95.control.nmi[[out]][[predictor]][,1]
  width.nmi[[out]][[predictor]] <- ci.95.nmi[[out]][[predictor]][,2] - ci.95.nmi[[out]][[predictor]][,1]
  
  q.bar.control.smi[[out]][[predictor]] <- apply(simplify2array(q.hat.control.smi[[out]][[predictor]]), 1, mean)
  q.bar.smi[[out]][[predictor]] <- apply(simplify2array(q.hat.smi[[out]][[predictor]]), 1, mean)
  var.u.control.smi[[out]][[predictor]] <- diag(apply(simplify2array(u.control.smi[[out]][[predictor]]), c(1:2), mean))
  var.u.smi[[out]][[predictor]] <- diag(apply(simplify2array(u.smi[[out]][[predictor]]), c(1:2), mean))
  e.b.control.smi[[out]][[predictor]] <- t(simplify2array(lapply(q.hat.control.smi[[out]][[predictor]], function(x) x - q.bar.control.smi[[out]][[predictor]])))
  e.b.smi[[out]][[predictor]] <- t(simplify2array(lapply(q.hat.smi[[out]][[predictor]], function(x) x - q.bar.smi[[out]][[predictor]])))
  var.b.control.smi[[out]][[predictor]] <- diag(t(e.b.control.smi[[out]][[predictor]])%*%e.b.control.smi[[out]][[predictor]] * m2.smi/(m1-1))
  var.b.smi[[out]][[predictor]] <- diag(t(e.b.smi[[out]][[predictor]])%*%e.b.smi[[out]][[predictor]] * m2.smi/(m1-1))
  var.w.control.smi[[out]][[predictor]] <- rep(0, length(q.bar.control.smi[[out]][[predictor]]))
  var.w.smi[[out]][[predictor]] <- rep(0, length(q.bar.smi[[out]][[predictor]]))
  var.t.control.smi[[out]][[predictor]] <- var.u.control.smi[[out]][[predictor]] + 1/m2.smi * (1 + 1/m1) * var.b.control.smi[[out]][[predictor]] + (1 - 1/m2.smi) * var.w.control.smi[[out]][[predictor]]
  var.t.smi[[out]][[predictor]] <- var.u.smi[[out]][[predictor]] + 1/m2.smi * (1 + 1/m1) * var.b.smi[[out]][[predictor]] + (1 - 1/m2.smi) * var.w.smi[[out]][[predictor]]
  df.control.smi[[out]][[predictor]] <- df.function(var.u.control.smi[[out]][[predictor]], var.b.control.smi[[out]][[predictor]], var.w.control.smi[[out]][[predictor]], var.t.control.smi[[out]][[predictor]], m1=m1, m2=m2.smi, 
                                              nrow(imp.data.nest.smi), length(q.bar.control.smi[[out]][[predictor]]), length(q.bar.control.smi[[out]][[predictor]]))
  df.smi[[out]][[predictor]] <- df.function(var.u.smi[[out]][[predictor]], var.b.smi[[out]][[predictor]], var.w.smi[[out]][[predictor]], var.t.smi[[out]][[predictor]], m1=m1, m2=m2.smi, 
                                        nrow(imp.data.nest.smi), length(q.bar.smi[[out]][[predictor]]), length(q.bar.smi[[out]][[predictor]]))
  q.t.control.smi[[out]][[predictor]] <- qt(0.975, df.control.smi[[out]][[predictor]])
  q.t.smi[[out]][[predictor]] <- qt(0.975, df.smi[[out]][[predictor]])
  t.value.control.smi[[out]][[predictor]] <- (q.bar.control.smi[[out]][[predictor]] - 0)/sqrt(var.t.control.smi[[out]][[predictor]])
  t.value.smi[[out]][[predictor]] <- (q.bar.smi[[out]][[predictor]] - 0)/sqrt(var.t.smi[[out]][[predictor]])
  p.value.control.smi[[out]][[predictor]] <- round(pt(t.value.control.smi[[out]][[predictor]], df.control.smi[[out]][[predictor]], lower.tail=FALSE), 3)
  p.value.smi[[out]][[predictor]] <- round(pt(t.value.smi[[out]][[predictor]], df.smi[[out]][[predictor]], lower.tail=FALSE), 3)
  ci.95.control.smi[[out]][[predictor]] <- round(cbind(q.bar.control.smi[[out]][[predictor]] - q.t.control.smi[[out]][[predictor]] * sqrt(var.t.control.smi[[out]][[predictor]]), 
                                                 q.bar.control.smi[[out]][[predictor]] + q.t.control.smi[[out]][[predictor]] * sqrt(var.t.control.smi[[out]][[predictor]])), 3)
  ci.95.smi[[out]][[predictor]] <- round(cbind(q.bar.smi[[out]][[predictor]] - q.t.smi[[out]][[predictor]] * sqrt(var.t.smi[[out]][[predictor]]), 
                                           q.bar.smi[[out]][[predictor]] + q.t.smi[[out]][[predictor]] * sqrt(var.t.smi[[out]][[predictor]])), 3)
  width.control.smi[[out]][[predictor]] <- ci.95.control.smi[[out]][[predictor]][,2] - ci.95.control.smi[[out]][[predictor]][,1]
  width.smi[[out]][[predictor]] <- ci.95.smi[[out]][[predictor]][,2] - ci.95.smi[[out]][[predictor]][,1]},
  warning=function(w)
    {message <- paste("Warning in predictor ", predictors.catchup[predictor], "\n", w)},
  error=function(e)
    {message <- paste("Error in predictor ", predictors.catchup[predictor], "\n", e)})}}

#### Save workspace ####
save.image("Workspaces/Workspace_POPS.Rdata")
