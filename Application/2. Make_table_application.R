#### POPS study ####
# Make sure folder "Application" is set as working directory

#### Basics ####
load("Workspaces/Workspace_POPS.Rdata")
require(mice)
require(xtable)
require(MASS)
require(extrafont)

loadfonts(device = "win")
par(family = "LM Roman 10")

#### Correlations ####
align.vector.corr <- c("lP{0.07\textwidth}P{0.07\textwidth}P{0.07\textwidth}
				P{0.07\textwidth}P{0.07\textwidth}P{0.07\textwidth}
				P{0.07\textwidth}P{0.07\textwidth}P{0.07\textwidth}")
addtorow.corr <- list()
addtorow.corr$pos <- list(0)
addtorow.corr$command <- c("& \\multicolumn{1}{l}{$\\text{Weight}_{1yr}$} & 
                   \\multicolumn{1}{l}{$\\text{Length}_{1yr}$} &  
                   \\multicolumn{1}{l}{$\\text{HC}_{1yr}$} & 
                   \\multicolumn{1}{l}{$\\text{WL}_{1yr}$} & 
                   \\multicolumn{1}{l}{Cognition} & 
                   \\multicolumn{1}{l}{HRQoL} & 
                   \\multicolumn{1}{l}{Int. beh.} & 
                   \\multicolumn{1}{l}{Ext. beh.} &
                   \\multicolumn{1}{l}{$\text{Length}_{19yr}$} \\\\ ")

print.xtable(xtable(upper.cor.pops, align.vector=align.vector.corr,
                    caption="Correlations between catch-up growth and outcomes 
             at age 19 among full responders in the POPS dataset.",
                    label="tab:cor_pops"), 
             include.rownames=TRUE, include.colnames=FALSE, add.to.row = addtorow.corr,
             caption.placement="top", NA.string=getOption("xtable.NA.string", " "))

#### Missingness ####
align.vector.miss <- c("lrrrrrrrrr")
addtorow.miss <- list()
addtorow.miss$pos <- list(0, 4)
addtorow.miss$command <- c("\\multicolumn{9}{l}{Unadjusted} \\\\
                      & \\multicolumn{1}{l}{$n_{\\text{total}$} & 
                      \\multicolumn{1}{l}{$n_{\\text{com}}$} &
                      \\multicolumn{1}{l}{\\% missing}&
                      \\multicolumn{1}{l}{Cognition}&
                      \\multicolumn{1}{l}{HRQol}&
                      \\multicolumn{1}{l}{Int. beh.}&
                      \\multicolumn{1}{l}{Ext. beh.}&
                      \\multicolumn{1}{l}{$\\text{Length}_{19yr}$} \\\\",
                           "\\multicolumn{9}{l}{Adjusted} \\\\
                      & \\multicolumn{1}{l}{$n_{\\text{total}$} & 
                      \\multicolumn{1}{l}{$n_{\\text{com}}$} &
                      \\multicolumn{1}{l}{\\% missing}&
                      \\multicolumn{1}{l}{Cognition}&
                      \\multicolumn{1}{l}{HRQol}&
                      \\multicolumn{1}{l}{Int. beh.}&
                      \\multicolumn{1}{l}{Ext. beh.}&
                      \\multicolumn{1}{l}{$\\text{Length}_{19yr}$} \\\\")

print(xtable(rbind(miss.df.unadj, miss.df.adj), align = align.vector.miss,
             caption="Missingness per data selection and dropout per analysis model in 
             the POPS data. Dropout is specified per outcome variable as a percentage 
             of total missingness.",
             label="tab:missingness"), 
      add.to.row=addtorow.miss,
      include.rownames = F, include.colnames = F,
      caption.placement="top", NA.string=getOption("xtable.NA.string", " "))


#### Results ####
tab.unadjusted <- tab.adjusted <- tab.combined <-
  relative.width <- relative.width.control <- vector("list", length(outcome))

# Restructure results 
for(out in 1:length(outcome)){
  q.bar.reimp.df <- rapply(q.bar.reimp[[out]], f=`[[`, ...=2, how="unlist")
  q.bar.nmi.df <- rapply(q.bar.nmi[[out]], f=`[[`, ...=2, how="unlist")
  q.bar.smi.df <- rapply(q.bar.smi[[out]], f=`[[`, ...=2, how="unlist")
  
  q.bar.control.reimp.df <- rapply(q.bar.control.reimp[[out]], f=`[[`, ...=2, how="unlist")
  q.bar.control.nmi.df <- rapply(q.bar.control.nmi[[out]], f=`[[`, ...=2, how="unlist")
  q.bar.control.smi.df <- rapply(q.bar.control.smi[[out]], f=`[[`, ...=2, how="unlist")
  
  
  ci.95.reimp.df <- t(simplify2array(lapply(ci.95.reimp[[out]],"[",2,)))
  ci.95.nmi.df <- t(simplify2array(lapply(ci.95.nmi[[out]],"[",2,)))
  ci.95.smi.df <- t(simplify2array(lapply(ci.95.smi[[out]],"[",2,)))
  
  ci.95.control.reimp.df <- t(simplify2array(lapply(ci.95.control.reimp[[out]],"[",2,)))
  ci.95.control.nmi.df <- t(simplify2array(lapply(ci.95.control.nmi[[out]],"[",2,)))
  ci.95.control.smi.df <- t(simplify2array(lapply(ci.95.control.smi[[out]],"[",2,)))
  
  width.reimp.df <- rapply(width.reimp[[out]], f=`[[`, ...=2, how="unlist")
  width.nmi.df <- rapply(width.nmi[[out]], f=`[[`, ...=2, how="unlist")
  width.smi.df <- rapply(width.smi[[out]], f=`[[`, ...=2, how="unlist")
  
  width.control.reimp.df <- rapply(width.control.reimp[[out]], f=`[[`, ...=2, how="unlist")
  width.control.nmi.df <- rapply(width.control.nmi[[out]], f=`[[`, ...=2, how="unlist")
  width.control.smi.df <- rapply(width.control.smi[[out]], f=`[[`, ...=2, how="unlist")
  
  #### Relative width ####
  relative.width[[out]] <- list(width.nmi.df/width.reimp.df, 
                                width.nmi.df/width.smi.df, 
                                width.smi.df/width.reimp.df)
  relative.width.control[[out]] <- list(width.control.nmi.df/width.control.reimp.df, 
                                        width.control.nmi.df/width.control.smi.df, 
                                        width.control.smi.df/width.control.reimp.df)
  #### Make table ####
  unadjusted.df <- cbind(formatC(q.bar.reimp.df, digits = 2, format = "f"),
                         formatC(ci.95.reimp.df, digits = 2, format = "f"),
                         formatC(width.reimp.df, digits = 2, format = "f"),
                         formatC(q.bar.nmi.df, digits = 2, format = "f"),
                         formatC(ci.95.nmi.df, digits = 2, format = "f"),
                         formatC(width.nmi.df, digits = 2, format = "f"),
                         formatC(q.bar.smi.df, digits = 2, format = "f"),
                         formatC(ci.95.smi.df, digits = 2, format = "f"),
                         formatC(width.smi.df, digits = 2, format = "f"))
  
  adjusted.df <- cbind(formatC(q.bar.control.reimp.df, digits = 2, format = "f"),
                       formatC(ci.95.control.reimp.df, digits = 2, format = "f"),
                       formatC(width.control.reimp.df, digits = 2, format = "f"),
                       formatC(q.bar.control.nmi.df, digits = 2, format = "f"),
                       formatC(ci.95.control.nmi.df, digits = 2, format = "f"),
                       formatC(width.control.nmi.df, digits = 2, format = "f"),
                       formatC(q.bar.control.smi.df, digits = 2, format = "f"),
                       formatC(ci.95.control.smi.df, digits = 2, format = "f"),
                       formatC(width.control.smi.df, digits = 2, format = "f"))
  unadjusted.df <- cbind(c("Weight", "Length", "Head circumference", "Weight-length"), unadjusted.df) 
  adjusted.df <- cbind(c("Weight", "Length", "Head circumference", "Weight-length"), adjusted.df) 
  
  combined.df <- rbind(unadjusted.df, adjusted.df)

  align.vector <- c("ll|cccc||cccc||cccc")
  addtorow.combined <- list()
  addtorow.combined$pos <- list()
  addtorow.combined$pos[[1]] <- 0
  addtorow.combined$pos[[2]] <- 4
  addtorow.combined$command <- c("\\multicolumn{13}{c}{Unadjusted} \\\\ \\hline
                                 Type catchup growth & \\multicolumn{4}{l}{Reimputation} & \\multicolumn{4}{l}{NMI} & \\multicolumn{4}{l}{SMI}\\\\ \\hline
                                 
                                 
                                 & B & \\multicolumn{2}{l}{95\\% CI} & Width & B & \\multicolumn{2}{l}{95\\% CI} & Width & 
                                 B & \\multicolumn{2}{l}{95\\% CI} & Width \\\\ \\hline",
                                 "\\hline \\multicolumn{13}{c}{Adjusted} \\\\ \\hline
                                 Type catchup growth & \\multicolumn{4}{l}{Reimputation} & \\multicolumn{4}{l}{NMI} & \\multicolumn{4}{l}{SMI}\\\\ \\hline
                                 
                                 
                                 & B & \\multicolumn{2}{l}{95\\% CI} & Width & B & \\multicolumn{2}{l}{95\\% CI} & Width & 
                                 B & \\multicolumn{2}{l}{95\\% CI} & Width \\\\ \\hline ")
  
  tab.combined[[out]] <- xtable(combined.df, align = align.vector)}


#### Presented in paper ####
  # Correlations
print.xtable(xtable(upper.cor.pops, align.vector=align.vector.corr,
                    caption="Correlations between catch-up growth and outcomes 
             at age 19 among full responders in the POPS dataset.",
                    label="tab:cor_pops"), 
             include.rownames=TRUE, include.colnames=FALSE, add.to.row = addtorow.corr,
             caption.placement="top", NA.string=getOption("xtable.NA.string", " "))
  # Missingness
print(xtable(rbind(miss.df.unadj, miss.df.adj), align = align.vector.miss,
             caption="Missingness per data selection and dropout per analysis model in 
             the POPS data. Dropout is specified per outcome variable as a percentage 
             of total missingness.",
             label="tab:missingness"), 
      add.to.row=addtorow.miss,
      include.rownames = F, include.colnames = F,
      caption.placement="top", NA.string=getOption("xtable.NA.string", " "))
  # Results
  print.xtable(tab.combined[[5]], add.to.row = addtorow.combined, include.rownames = F, 
        include.colnames = F)

  ## Relative width  
  # (1: NMI vs. RMI; 2: NMI vs. AMI; 3: NMI vs. SMI)
  # for Weight, Length, Head circumference, Weight-length
  # Unadjusted
(relative.width.unadj <- relative.width[[5]])
  # Adjusted 
(relative.width.adj <- relative.width.control[[5]])