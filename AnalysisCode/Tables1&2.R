###########################################
### Spatial Clustering - Tables 1-2   #####
###########################################


library(reshape2)
library(SpatialClusteringPkg)

source("./R/outbreak_prop_functs.R")
source("./FinalSize/Source/finalsize_utilities.R")


write.excel <- function(x,row.names=TRUE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}

# multi.sir <- function(ntrials=10, n=1000,...){
#   res <- foreach(i=1:ntrials, .errorhandling="pass", .combine=cbind, .multicombine=TRUE) %do% {
#     simulateFinalSizes(n, ...)
#   }
#   return(res)
# }





#*********************************************************************
#*********************************************************************
#   TABLE 1                                  ##### *******************
#   - Reff analyitic estimates using multiple g(x) distributions
#*********************************************************************
VE=1

#Read et al. 2014 Contact Distrib #

upper.lim=400
vacc <- c(.8,.85,.90,.91,.92,.93,.94,.95)
thetas <- c(0, .25, .5, 1) # Theta = 0 (no clustering), 0.25 (low), 0.5, (medium), and 1.0 (high) 


get.tau.params

A <- sapply(X=vacc, FUN=function(Y=X) sapply(X=thetas, FUN=calc.R.cf, 
                                             v=Y, VE=VE, cbeta=15, lambda=.5, psi=1,
                                             alpha=g_params$china$shape, beta=1/g_params$china$scale))
colnames(A) <- vacc; row.names(A) <- thetas+1
A


#Zambia Contact Distrib #
B <- sapply(X=vacc, FUN=function(Y=X) sapply(X=thetas, FUN=calc.R.cf, 
                                             v=Y, VE=VE, cbeta=15, lambda=.5, psi=1,
                                             alpha=g_params$zambia$shape, beta=1/g_params$zambia$scale))
colnames(B) <- vacc; row.names(B) <- thetas+1
B


#Noulas et al. 2012 Contact Distrib #
C <- sapply(X=vacc, FUN=function(Y=X) sapply(X=thetas, FUN=calc.R.cf, 
                                             v=Y, VE=VE, cbeta=15, lambda=.5, psi=1,
                                             alpha=g_params$usa$shape, beta=1/g_params$usa$scale))
colnames(C) <- vacc; row.names(C) <- thetas+1
C

tau_levels <- c("None", "Low", "Medium", "High")
vacc_levels <- c("80%","","","","85%","","","","90%","","","","91%","","","","92%","","","","93%","","","","94%","","","","95%","","","")
tab1 <- cbind(vacc_levels, rep(tau_levels,length(vacc)), round(as.vector(A),2), round(as.vector(B),2), round(as.vector(C),2))
colnames(tab1) <- c("Vaccination Coverage", "Clustering", "Read et al. 2014", "Zambia", "Noulas et al. 2012")
tab1
#write.excel(tab1,row.names=FALSE,col.names=TRUE)


# SAVE R ESTIMATES TO FILE
R_ests <- tab1
R_ests[,1] <- sort(rep(vacc, length(thetas)))
R_ests[,2] <- rep(thetas+1, length(vacc))
R_ests <- matrix(as.numeric(R_ests), ncol=4, nrow=length(vacc)*length(thetas))
colnames(R_ests) <- c("vacc", "tau", "Read", "Noulas")
save(R_ests, file="./R_estimates.RData")



# Easy to view version
R_tab <- data.frame(gA=A[,'0.95'], gB=B[,'0.95'], gC=C[,'0.95'])
R_tab









#*********************************************************************
#   - Vc analyitic estimates using multiple g(x) distributions

VE = 1
# Read et al. 2014 Contact Distrib #

# Theta = 0 (no clustering), 0.25 (low), 0.5, (medium), and 1.0 (high) 
VcA <- sapply(X=c(0, .25, .5, 1), FUN=calc.Vc.cf, 
              VE=VE, R0=15, lambda=.5, psi=1,
              alpha=g_params$china$shape, beta=1/g_params$china$scale)
VcA

# Noulas et al. 2012 Contact Distrib #
# Theta = 0 (no clustering), 0.25 (low), 0.5, (medium), and 1.0 (high) 
VcB <- sapply(X=c(0, .25, .5, 1), FUN=calc.Vc.cf, 
              VE=VE, R0=15, lambda=.5, psi=1,
              alpha=g_params$zambia$shape, beta=1/g_params$zambia$scale)
VcB

VcC <- sapply(X=c(0, .25, .5, 1), FUN=calc.Vc.cf, 
              VE=VE, R0=15, lambda=.5, psi=1,
              alpha=g_params$usa$shape, beta=1/g_params$usa$scale)
VcC


tau_levels <- c("None", "Low", "Medium", "High")
tab2 <- cbind(tau_levels, round(VcA,4), round(VcB,4), round(VcC,4))
colnames(tab2) <- c("Clustering", "Read et al. 2014", "Zambia", "Noulas et al. 2012")
tab2

#write.excel(tab2,row.names=FALSE,col.names=TRUE)

Vc_tab <- data.frame(gA=VcA, gB=VcB, gC=VcC)
Vc_tab




R_tab
Vc_tab








#*********************************************************************
#*********************************************************************
#   TABLE 2                                  ##### *******************
#   - Deterministic, SIR, and Spatial Estimates of Outbreak Probability
#*********************************************************************

# SAVE R ESTIMATES TO FILE
R_ests <- tab1
R_ests[,1] <- sort(rep(vacc, length(thetas)))
R_ests[,2] <- rep(thetas+1, length(vacc))
R_ests <- matrix(as.numeric(R_ests), ncol=5, nrow=length(vacc)*length(thetas))
colnames(R_ests) <- c("vacc", "tau", "Read", "zambia", "Noulas")
save(R_ests, file="./Results/R_estimates.RData")

load(file="./Results/R_estimates.RData")





# Sim - 80%    ####
Sim80pt <- outbreak.prop.sims(tau=2, vacc.pt=.80, spatial=TRUE, plot.hist=TRUE, folder="step.size_1", cutoff=.05,
                              popfolder="30x25km", i=1, gridsize=30, pop.dens=25, upper.lim=500, lambda=0.4); Sim80pt
res <- sapply(c(1.25, 1.5, 2.0), outbreak.prop.sims, vacc.pt=.80, spatial=TRUE, plot.hist=TRUE, folder="step.size_1", cutoff=.05,
              popfolder="30x25km", i=1, gridsize=30, pop.dens=25, upper.lim=500, lambda=0.4); t(res)
write.excel(t(res)[,-c(1,2,7,8)], row.names = F, col.names = F)

# Sim - 85%    ####
Sim85pt <- outbreak.prop.sims(tau=1.5, vacc.pt=.85, spatial=TRUE, plot.hist=TRUE, folder="step.size_1", cutoff=.05,
                              popfolder="30x25km", i=1, gridsize=30, pop.dens=25, upper.lim=500, lambda=0.4); Sim85pt
res <- sapply(c(1.25, 1.5, 2.0), outbreak.prop.sims, vacc.pt=.90, spatial=TRUE, plot.hist=F, folder="step.size_1", cutoff=.05,
              popfolder="30x25km", i=1, gridsize=30, pop.dens=25, upper.lim=500, lambda=0.4); t(res)
write.excel(t(res)[,-c(1,2,7,8)], row.names = F, col.names = F)

# Sim - 90%    ####
Sim90pt <- outbreak.prop.sims(tau=1.5, vacc.pt=.90, spatial=TRUE, plot.hist=TRUE, folder="step.size_1", cutoff=.05,
                              popfolder="30x25km", i=1, gridsize=30, pop.dens=25, upper.lim=500, lambda=0.4); Sim90pt
res <- sapply(c(1.25, 1.5, 2.0), outbreak.prop.sims, vacc.pt=.90, spatial=TRUE, plot.hist=F, folder="step.size_1", cutoff=.05,
              popfolder="30x25km", i=1, gridsize=30, pop.dens=25, upper.lim=500, lambda=0.4); t(res)
write.excel(t(res)[,-c(1,2,7,8)], row.names = F, col.names = F)

# Sim - 91%    ####
res <- sapply(c(1.25, 1.5, 2.0), outbreak.prop.sims, vacc.pt=.91, spatial=TRUE, plot.hist=F, folder="step.size_1", cutoff=.05,
              popfolder="30x25km", i=1, gridsize=30, pop.dens=25, upper.lim=500, lambda=0.4); t(res)
write.excel(t(res)[,-c(1,2,7,8)], row.names = F, col.names = F)

# Sim - 92%    ####
res <- sapply(c(1.25, 1.5, 2.0), outbreak.prop.sims, vacc.pt=.92, spatial=TRUE, plot.hist=TRUE, folder="step.size_1", cutoff=.05,
              popfolder="30x25km", i=1, gridsize=30, pop.dens=25, upper.lim=500, lambda=0.4); t(res)
write.excel(t(res)[,-c(1,2,7,8)], row.names = F, col.names = F)

# Sim - 93%    ####
res <- sapply(c(1.25, 1.5, 2.0), outbreak.prop.sims, vacc.pt=.93, spatial=TRUE, plot.hist=F, folder="step.size_1", cutoff=.05,
              popfolder="30x25km", i=1, gridsize=30, pop.dens=25, upper.lim=500, lambda=0.4); t(res)
write.excel(t(res)[,-c(1,2,7,8)], row.names = F, col.names = F)

# Sim - 94%    ####
res <- sapply(c(1.25, 1.5, 2.0), outbreak.prop.sims, vacc.pt=.94, spatial=TRUE, plot.hist=TRUE, folder="step.size_1", cutoff=.05,
              popfolder="30x25km", i=1, gridsize=30, pop.dens=25, upper.lim=500, lambda=0.4); t(res)
write.excel(t(res)[,-c(1,2,7,8)], row.names = F, col.names = F)

# Sim - 95%    ####
res <- sapply(c(1.25, 1.5, 2.0), outbreak.prop.sims, vacc.pt=.95, spatial=TRUE, plot.hist=TRUE, folder="step.size_1", cutoff=.05,
              popfolder="30x25km", i=1, gridsize=30, pop.dens=25, upper.lim=500, lambda=0.4); t(res)
write.excel(t(res)[,-c(1,2,7,8)], row.names = F, col.names = F)





#--- SIR Simulations ---#
load(file="./Results/R_estimates.RData") # Loads "reff_ests"

# SIR - 80%    ####
SIR80pt <- outbreak.prop.SIR(vacc.pt=.80, step.size=.1, plot.hist=TRUE, folder="SIR", cutoff=.05, upper.lim=500, lambda=0.4); SIR80pt
write.excel(SIR80pt, row.names = F, col.names = F)

# SIR - 85%    ####
SIR85pt <- outbreak.prop.SIR(vacc.pt=.85, step.size=.1, plot.hist=TRUE, folder="SIR", cutoff=.05, upper.lim=500, lambda=0.4); SIR85pt
write.excel(SIR85pt, row.names = F, col.names = F)

# SIR - 90%    ####
SIR90pt <- outbreak.prop.SIR(vacc.pt=.90, step.size=.1, plot.hist=TRUE, folder="SIR", cutoff=.05, upper.lim=500, lambda=0.4); SIR90pt
write.excel(SIR90pt, row.names = F, col.names = F)

# SIR - 91%    ####
SIR91pt <- outbreak.prop.SIR(vacc.pt=.91, step.size=.1, plot.hist=TRUE, folder="SIR", cutoff=.05, upper.lim=500, lambda=0.4); SIR91pt
write.excel(SIR91pt, row.names = F, col.names = F)

# SIR - 92%    ####
SIR92pt <- outbreak.prop.SIR(vacc.pt=.92, step.size=.1, plot.hist=TRUE, folder="SIR", cutoff=.05, upper.lim=500, lambda=0.4); SIR92pt
write.excel(SIR92pt, row.names = F, col.names = F)

# SIR - 93%    ####
SIR93pt <- outbreak.prop.SIR(vacc.pt=.93, step.size=.1, plot.hist=TRUE, folder="SIR", cutoff=.05, upper.lim=500, lambda=0.4); SIR93pt
write.excel(SIR93pt, row.names = F, col.names = F)

# SIR - 94%    ####
SIR94pt <- outbreak.prop.SIR(vacc.pt=.94, step.size=.1, plot.hist=TRUE, folder="SIR", cutoff=.05, upper.lim=500, lambda=0.4); SIR94pt
write.excel(SIR94pt, row.names = F, col.names = F)

# SIR - 95%    ####
SIR95pt <- outbreak.prop.SIR(vacc.pt=.95, step.size=.1, plot.hist=TRUE, folder="SIR", cutoff=.05, upper.lim=500, lambda=0.4); SIR95pt
write.excel(SIR95pt, row.names = F, col.names = F)



#### Correlations ####
corr.funct <- function(vacc.pt, cutoff,...){
  combdat <- outbreak.prop.all(vacc.pt=vacc.pt, cutoff=cutoff, folder="step.size_1", popfolder="30x25km", step.size=1, step.size.sir=.1, sirfolder="SIR")
  tmp <- cbind(sir=combdat[1:4,"mean"], sim=combdat[5:8, "mean"])
  return(c(corr=cor(tmp[,1], tmp[,2]), chisq=chisq.test(tmp)$p.value))
}

stat.tests <- function(vacc.pt,cutoff,...){
  combdat <- outbreak.prop.all(vacc.pt=vacc.pt, cutoff=cutoff,folder="step.size_1", popfolder="30x25km", step.size=1, step.size.sir=.1, sirfolder="SIR")
  tmp <- cbind(sir=combdat[1:4,"mean"], sim=combdat[5:8, "mean"])
  tmp2 <- cbind(type=c(0,0,0,0,1,1,1,1), mean=combdat$mean)
  MannWhit <- wilcox.test(tmp2[,2]~tmp2[,1])
  
  return(c(corr=cor(tmp[,1], tmp[,2]), chisq=chisq.test(tmp)$p.value))
}

proport.funct <- function(vacc.pt, cutoff,...){
  combdat <- outbreak.prop.all(vacc.pt=vacc.pt, cutoff=cutoff, folder="step.size_1", popfolder="30x25km", step.size=1, step.size.sir=.1, sirfolder="SIR")
  tmp <- cbind(sir=combdat[1:4,"mean"], sim=combdat[5:8, "mean"])
  return(mean(tmp[,2] / tmp[,1]))
}

cutoff=.05
all95pt <- outbreak.prop.all(vacc.pt=.95, cutoff=cutoff, folder="step.size_1", popfolder="30x25km", step.size=1, step.size.sir=.1, sirfolder="SIR")
all94pt <- outbreak.prop.all(vacc.pt=.94, cutoff=cutoff, folder="step.size_1", popfolder="30x25km", step.size=1, step.size.sir=.1, sirfolder="SIR")
all93pt <- outbreak.prop.all(vacc.pt=.93, cutoff=cutoff, folder="step.size_1", popfolder="30x25km", step.size=1, step.size.sir=.1, sirfolder="SIR")
all92pt <- outbreak.prop.all(vacc.pt=.92, cutoff=cutoff, folder="step.size_1", popfolder="30x25km", step.size=1, step.size.sir=.1, sirfolder="SIR")
all91pt <- outbreak.prop.all(vacc.pt=.91, cutoff=cutoff, folder="step.size_1", popfolder="30x25km", step.size=1, step.size.sir=.1, sirfolder="SIR")
all90pt <- outbreak.prop.all(vacc.pt=.90, cutoff=cutoff, folder="step.size_1", popfolder="30x25km", step.size=1, step.size.sir=.1, sirfolder="SIR")
all85pt <- outbreak.prop.all(vacc.pt=.85, cutoff=cutoff, folder="step.size_1", popfolder="30x25km", step.size=1, step.size.sir=.1, sirfolder="SIR")
all80pt <- outbreak.prop.all(vacc.pt=.80, cutoff=cutoff, folder="step.size_1", popfolder="30x25km", step.size=1, step.size.sir=.1, sirfolder="SIR")

results.all <- rbind(all95pt, all94pt, all93pt, all92pt, all91pt, all90pt, all80pt)

# Plot of Divergence of Outbreak Probability #
probplot1 <- ggplot(results.all, aes(x = vacc, y = mean, color=type, linetype=as.factor(tau))) + 
  geom_line(size=1.25) + scale_color_manual(values=c("blue", "red")) +
  scale_linetype_manual(values=c("solid", "solid", 'solid', 'solid'))
probplot1 <- probplot1 + ylab("E[Pr(outbreak)]") + xlab("Effective Vaccination Coverage")
probplot1 <- probplot1 + theme_bw() + theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
print(probplot1)



# Correlation for Paper #
corrdat <- cbind(results.all[which(results.all$type=="sir"),"mean"], results.all[which(results.all$type=="sims"),"mean"])
cor(corrdat[,1], corrdat[,2], method="pearson")
cor(corrdat[,1], corrdat[,2], method="spearman")

chisq.test(y=corrdat[,1], x=corrdat[,2])




lapply(seq(.91, .94, .01), corr.funct, cutoff=.05)
lapply(seq(.91, .95, .01), proport.funct, cutoff=.05)

corr.funct(vacc.pt=.95,cutoff=.05)
corr.funct(vacc.pt=.93,cutoff=.05)

#Correlations
vacc <- c(.8,.90,.91,.92,.93,.94,.95)
lapply(vacc, corr.funct, cutoff=.05)



proport.funct(vacc.pt=.91,cutoff=.05)

tmp <- aggregate(alldat, by=list(as.factor(alldat$type), as.factor(alldat$reff)), FUN="mean")

model1 <- glm(prop ~ type + reff, data=alldat, family="binomial")




all.95pt.check <- outbreak.check.all(vacc.pt=.95, cutoff=.05)
model1 <- glm(outbreak ~ type*reff, data=all.95pt.check, family="binomial")
summary(model1)
exp(coef(model1))

all.93pt.check <- outbreak.check.all(vacc.pt=.93, cutoff=.05)
model1 <- glm(outbreak ~ type*reff, data=all.93pt.check, family="binomial")
summary(model1)
exp(coef(model1))

all.92pt.check <- outbreak.check.all(vacc.pt=.92, cutoff=.10)
model1 <- glm(outbreak ~ type*reff, data=all.92pt.check, family="binomial")
summary(model1)
exp(coef(model1))


all.92pt.check <- outbreak.check.all(vacc.pt=.92, cutoff=.10)
reff92 <- unique(all.92pt.check$reff)
model1 <- glm(outbreak ~ type, data=all.92pt.check[which(all.92pt.check$reff==reff92[4]),], family="binomial")
summary(model1)
exp(coef(model1))














