#### POPS study ####
rm(list=ls())
#### Basics ####
require(mice)
require(xtable)
require(MASS)

set.seed(13579)

m1 <- 5


source("Df.R")


#### CATCHUP GROWTH ####
# Select necessary variables
pops.totaal <- read.csv("POPS_totaal_numeric.csv", header=TRUE)
colnames(pops.totaal)[1] <- "laatste"
outcomes <- c("btchol", "btc", "bhchol", "aa7", "bmistatus_meancat2rev",
              "iq", "LHS_POPS19x", "MAU_POPS19xcat4",   "Internal_Geel", "External_Geel", 
              "bmistatus_meancat3", "l19_sd")
predict.weight <- c("G_SD_GA1yr", "G_SD_GA0yr")
predict.length <- c("L_SD_GA1yr", "L_SD_GA0yr")
predictors.length <- c("L_SD_GA1yr", "L_SD_GA0yr", "l5_sd", "l10_sd", "l14_sd")
predictors.weight <- c("G_SD_GA1yr", "G_SD_GA0yr", "g5_sd", "g10_sd", "g14_sd", "g19_sd")
predictors.control <- c("mag", "mag2", "rok", "var148cat2", "sestot", "thsds", "par", "zwduur", "zwduur2",
                        "geslr", "ras")
predictors.hc <- c("H_SD_GA1yr", "H_SD_GA0yr")
predictors.wl <- c("GL_SD_GA1yr", "GL_SD_GA0yr")
predictors.catchup <- c(predictors.weight, predictors.length, predictors.hc, predictors.wl)
predictors <- c(predictors.weight, predictors.length, predictors.hc, predictors.wl, 
                predictors.control)
predictor.names <- c("Weight", "Length", "Head circumference", "Weight-length")
selection.vars <- c("sgalen", "sgagew", "sgaho", "sgagewlen")
vars <- c("patid", outcomes, predictors)
pops.data <- pops.totaal[,vars] 
pops.data <- apply(pops.data, 2, as.numeric)

selection.data <- pops.totaal[,selection.vars]
selection.data.long <- selection.data[rep(1:nrow(selection.data), times = m1),]

#### Data ####
# Available data at each wave
pops.0 <- pops.data[,c("patid", "L_SD_GA0yr", "G_SD_GA0yr", "GL_SD_GA0yr","H_SD_GA0yr", predictors.control)]
pops.1 <- cbind(pops.0, pops.data[,c("L_SD_GA1yr", "G_SD_GA1yr", "GL_SD_GA1yr","H_SD_GA1yr")])
pops.5 <- cbind(pops.1, pops.data[,c("l5_sd", "g5_sd")])
pops.10 <- cbind(pops.5, pops.data[,c("l10_sd", "g10_sd")])
pops.14 <- cbind(pops.10, pops.data[,c("l14_sd", "g14_sd")])
pops.19 <- cbind(pops.14, pops.data[,c("g19_sd")])
pops.outcomes <- pops.data[,outcomes]
pops.19 <- cbind(pops.19, pops.outcomes)

#### Pred ####
# Create predictor matrices for each age, based on available data at that wave
pred.0 <- matrix(rep(1, ncol(pops.0)^2), nrow=ncol(pops.0), ncol=ncol(pops.0), byrow=TRUE)
pred.5 <- matrix(rep(1, ncol(pops.5)^2), nrow=ncol(pops.5), ncol=ncol(pops.5), byrow=TRUE)
pred.10 <- matrix(rep(1, ncol(pops.10)^2), nrow=ncol(pops.10), ncol=ncol(pops.10), byrow=TRUE)
pred.14 <- matrix(rep(1, ncol(pops.14)^2), nrow=ncol(pops.14), ncol=ncol(pops.14), byrow=TRUE)
pred.19 <- matrix(rep(1, ncol(pops.19)^2), nrow=ncol(pops.19), ncol=ncol(pops.19), byrow=TRUE)
diag(pred.0) <- diag(pred.5) <- diag(pred.10) <- diag(pred.14) <- diag(pred.19) <- 0
colnames(pred.0) <- rownames(pred.0) <- colnames(pops.0)
colnames(pred.5) <- rownames(pred.5) <- colnames(pops.5)
colnames(pred.10) <- rownames(pred.10) <- colnames(pops.10)
colnames(pred.14) <- rownames(pred.14) <- colnames(pops.14)
colnames(pred.19) <- rownames(pred.19) <- colnames(pops.19)

#### Reimputation ####
reimp.0 <- mice(pops.0, predictorMatrix=pred.0, print=FALSE)
reimp.5 <- mice(pops.5, predictorMatrix=pred.5, print=FALSE)
reimp.10 <- mice(pops.10, predictorMatrix=pred.10, print=FALSE)
reimp.14 <- mice(pops.14, predictorMatrix=pred.14, print=FALSE)
reimp.19 <- mice(pops.19, predictorMatrix=pred.19, print=FALSE)

#### Weight ####
pops.weight.0.reimp.data <- cbind(complete(reimp.0, action="long", include=FALSE), selection.data.long)
pops.weight.5.reimp.data <- cbind(complete(reimp.5, action="long", include=FALSE), selection.data.long)
pops.weight.10.reimp.data <- cbind(complete(reimp.10, action="long", include=FALSE), selection.data.long)
pops.weight.14.reimp.data <- cbind(complete(reimp.14, action="long", include=FALSE), selection.data.long)
pops.weight.19.reimp.data <- cbind(complete(reimp.19, action="long", include=FALSE), selection.data.long)
pops.weight.0.reimp.data[,selection.vars] <- apply(pops.weight.0.reimp.data[,selection.vars], 2, as.numeric)
pops.weight.5.reimp.data[,selection.vars] <- apply(pops.weight.5.reimp.data[,selection.vars], 2, as.numeric)
pops.weight.10.reimp.data[,selection.vars] <- apply(pops.weight.10.reimp.data[,selection.vars], 2, as.numeric)
pops.weight.14.reimp.data[,selection.vars] <- apply(pops.weight.14.reimp.data[,selection.vars], 2, as.numeric)
pops.weight.19.reimp.data[,selection.vars] <- apply(pops.weight.19.reimp.data[,selection.vars], 2, as.numeric)

pops.weight.data <- list(pops.weight.0.reimp.data, pops.weight.5.reimp.data, pops.weight.10.reimp.data, 
                         pops.weight.14.reimp.data, pops.weight.19.reimp.data)
#### Weight: Analyses ####
# Predict length at age 5, 10, 14 and weight at age 5, 10, 14 from catch-up growth in weight
outcomes.lm <- c("l5_sd", "l10_sd", "l14_sd", "g5_sd", "g10_sd", "g14_sd") #Length at age 5, 10, 14; Weight at age 5, 10, 14
formulas.weight.control <- formulas.weight <-  vector("list", length(outcomes.lm))
# Create formulas
for(out in 1:length(outcomes.lm)){
  formulas.weight.control[[out]] <- as.formula(paste(outcomes.lm[out], paste(c(predict.weight, predictors.control), collapse=" + "), sep=" ~ "))
  formulas.weight[[out]] <- as.formula(paste(outcomes.lm[out], paste(predict.weight, collapse=" + "), sep=" ~ "))}

# Fit models and pool results
summary.weight.control.reimp <- summary.weight.reimp <- rep(list(list()), length(pops.weight.data)-1)
fit.control.reimp <- fit.reimp <- vector("list", m1)
for(wave in 2:length(pops.weight.data)){
  
  wave.data <- pops.weight.data[[wave]]
  for(out in 1:length(outcomes.lm)){
    tryCatch({
      for(imp.m1 in 1:5){
        imp.data <- wave.data[wave.data$.imp == imp.m1,]
        imp.data.reimp <- imp.data[imp.data[,"sgagew"] == 1 &
                                     !is.na(imp.data[,"sgagew"]) ,]
        fit.control.reimp[[imp.m1]] <- lm(formulas.weight.control[[out]], data = imp.data.reimp)
        fit.reimp[[imp.m1]] <- lm(formulas.weight[[out]], data = imp.data.reimp)}
      fit.control.reimp.mira <- as.mira(fit.control.reimp)
      fit.reimp.mira <- as.mira(fit.reimp)
      pool.control.reimp <- pool(fit.control.reimp.mira)
      pool.reimp <- pool(fit.reimp.mira)
      summary.weight.control.reimp[[wave-1]][[out]] <- summary(pool.control.reimp)
      summary.weight.reimp[[wave-1]][[out]] <- summary(pool.reimp)
    },
    warning=function(w)
    {message <- paste("Warning in outcome", outcomes.lm[out], "\n", w)
    print(message)},
    error=function(e)
    {message <- paste("Error in outcome", outcomes.lm[out], "\n", e)
    print(message)
    #summary.weight.control.reimp[[wave-1]][[out]] <- 
    return(matrix(1, nrow=14, ncol=9))
    #summary.weight.reimp[[wave-1]][[out]] <- matrix(1, nrow=3, ncol=9)
    })
  }}

#### Weight: Results ####
weight.5.1.df <- round(t(sapply(summary.weight.reimp, function(x) x[[1]][2,])), 3)
weight.5.2.df <- round(t(sapply(summary.weight.reimp, function(x) x[[4]][2,])), 3)
weight.10.1.df <-round(rbind(summary.weight.reimp[[2]][[2]][2,],
                             summary.weight.reimp[[3]][[2]][2,],
                             summary.weight.reimp[[4]][[2]][2,]), 3)
weight.10.2.df <-round(rbind(summary.weight.reimp[[2]][[5]][2,],
                             summary.weight.reimp[[3]][[5]][2,],
                             summary.weight.reimp[[4]][[5]][2,]), 3)
weight.14.1.df <- round(rbind(summary.weight.reimp[[3]][[3]][2,],
                              summary.weight.reimp[[4]][[3]][2,]), 3)
weight.14.2.df <- round(rbind(summary.weight.reimp[[3]][[6]][2,],
                              summary.weight.reimp[[4]][[6]][2,]), 3)

weight.control.5.1.df <- round(t(sapply(summary.weight.control.reimp, function(x) x[[1]][2,])), 3)
weight.control.5.2.df <- round(t(sapply(summary.weight.control.reimp, function(x) x[[4]][2,])), 3)
weight.control.10.1.df <-round(rbind(summary.weight.control.reimp[[2]][[2]][2,],
                                     summary.weight.control.reimp[[3]][[2]][2,],
                                     summary.weight.control.reimp[[4]][[2]][2,]), 3)
weight.control.10.2.df <-round(rbind(summary.weight.control.reimp[[2]][[5]][2,],
                                     summary.weight.control.reimp[[3]][[5]][2,],
                                     summary.weight.control.reimp[[4]][[5]][2,]), 3)
weight.control.14.1.df <- round(rbind(summary.weight.control.reimp[[3]][[3]][2,],
                                      summary.weight.control.reimp[[4]][[3]][2,]), 3)
weight.control.14.2.df <- round(rbind(summary.weight.control.reimp[[3]][[6]][2,],
                                      summary.weight.control.reimp[[4]][[6]][2,]), 3)

weight.5.1.df
weight.5.2.df
weight.10.1.df
weight.10.2.df
weight.14.1.df
weight.14.2.df

weight.control.5.1.df
weight.control.5.2.df
weight.control.10.1.df
weight.control.10.2.df
weight.control.14.1.df # Differs from 14.2 (data age 19)
weight.control.14.2.df # Differs from 14.1 (data age 14)




#### Presented in paper ####
data.age.14 <- weight.control.14.1.df[,c("est", "lo 95", "hi 95", "Pr(>|t|)")]
data.age.19 <- weight.control.14.2.df[,c("est", "lo 95", "hi 95", "Pr(>|t|)")]
length.by.weight <- rbind(data.age.14, data.age.19)
rownames(length.by.weight) <- c("Length - data age 14", "Weight - data age 14", "Length - data age 19", "Weight - data age 19")





