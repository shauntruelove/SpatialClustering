##' Figure 2
##'
##' Part 1:
##'
##' Part 2: Generate a heatmap examining the relationship between phi (the clustering factor) and
##'   Ic* (the clustering adjusted critical immunity coefficient) across R0 values.
##'
##'





### Paper1_Figure2 ######

# source("./R/clust_estimates_source.R")

library(SpatialClusteringPkg)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(reshape2)

# library(scales)
# require(graphics); require(grDevices); library(gplots); library(lattice)

## load Contact Parameters ####
# load(file="./data/g_params.rda")




#******************************************************************
#******************************************************************
####                     FIGURE 2A                             ####
#   - R0 vs Vc for all levels of clustering where R=1             #
#******************************************************************
#******************************************************************

VE=1
upper.lim=2000; lower.lim=0
vacc <- seq(.5, 1, .001)
thetas <- c(0, .25, .5, 1) # Theta = 0 (no clustering), 0.25 (low), 0.5, (medium), and 1.0 (high)
R0 <- seq(2, 25, .01)

# No longer need this. Closed form (below) is better)
# A <- sapply(X=R0, FUN=function(Y=X) sapply(X=thetas, FUN=calc.Vc.analytic, subdivisions=5000L, lower=lower.lim, by.x=0.001,
#                                            R0=Y, VE=VE, upper.lim=upper.lim, lambda=.5,
#                                            gx.funct.name="dgamma", shape=china.params$shape, scale=china.params$scale))
# colnames(A) <- R0; row.names(A) <- thetas+1
# A

A <- sapply(X=R0, FUN=function(Y=X) sapply(X=thetas, FUN=calc.Vc.cf,
                                           R0=Y, VE=VE, psi=1, lambda=.5,
                                           alpha=g_params$china$shape, beta=1/g_params$china$scale))
colnames(A) <- R0; row.names(A) <- thetas+1

Adata <- t(A)
Adata <- data.frame(cbind(row.names(Adata), Adata),  stringsAsFactors=F)
colnames(Adata) <- c("R0", "None", "Low", "Medium", "High")
Adata <- melt(Adata, id="R0")
colnames(Adata) <- c("R0", "Clustering", "Vc")
Adata$R0 <- as.numeric(as.character(Adata$R0))
Adata$Vc <- as.numeric(as.character(Adata$Vc))

R0plot1 <- ggplot(Adata, aes(x = R0, y = Vc, color=Clustering)) + geom_line(size=1.25) +
    scale_color_manual(values=c("blue", "gold1", "darkorange1", "red3")) +
    scale_x_continuous(trans="log", breaks=c(2, 5, 10, 15, 20), limits=c(2, 22))
R0plot1 <- R0plot1 + xlab(bquote('R'[0])) + ylab(bquote("Critical Vaccination Threshold (V"[c]~')'))
R0plot1 <- R0plot1 + theme_bw() + theme(axis.text=element_text(size=12),
                                        axis.title=element_text(size=14,face="bold"))
plot(R0plot1)


#*********************************************************************
# Zambia contact kernel  ####

B <- sapply(X=R0, FUN=function(Y=X) sapply(X=thetas, FUN=calc.Vc.cf,
                                           R0=Y, VE=VE, psi=1, lambda=.5,
                                           alpha=g_params$zambia$shape, beta=1/g_params$zambia$scale))
colnames(B) <- R0; row.names(B) <- thetas+1

Bdata <- t(B)
Bdata <- data.frame(cbind(row.names(Bdata), Bdata),  stringsAsFactors=F)
colnames(Bdata) <- c("R0", "None", "Low", "Medium", "High")
Bdata <- melt(Bdata, id="R0")
colnames(Bdata) <- c("R0", "Clustering", "Vc")
Bdata$R0 <- as.numeric(as.character(Bdata$R0))
Bdata$Vc <- as.numeric(as.character(Bdata$Vc))


# ------------------------------------------------------------------
# Noulas contact kernel  ####

C <- sapply(X=R0, FUN=function(Y=X) sapply(X=thetas, FUN=calc.Vc.cf,
                                           R0=Y, VE=VE, psi=1, lambda=.5,
                                           alpha=g_params$usa$shape, beta=1/g_params$usa$scale))
colnames(C) <- R0; row.names(C) <- thetas+1

Cdata <- t(C)
Cdata <- data.frame(cbind(row.names(Cdata), Cdata),  stringsAsFactors=F)
colnames(Cdata) <- c("R0", "None", "Low", "Medium", "High")
Cdata <- melt(Cdata, id="R0")
colnames(Cdata) <- c("R0", "Clustering", "Vc")
Cdata$R0 <- as.numeric(as.character(Cdata$R0))
Cdata$Vc <- as.numeric(as.character(Cdata$Vc))



######################################
# Combined Figure 3 - 2 plots - Log Scale X & Y  ####

combdata <- rbind(cbind(gx=rep("A", nrow(Adata)), Adata), cbind(gx=rep("B", nrow(Bdata)), Bdata),  cbind(gx=rep("C", nrow(Cdata)), Cdata))
combdata$gx <- as.factor(combdata$gx)
levels(combdata$gx) <- list('italic(g[A](r))'='A',
                            'italic(g[B](r))'='B',
                            'italic(g[C](r))'='C')
combdata$Clustering <- factor(combdata$Clustering, levels=rev(levels(combdata$Clustering)))

myPalette <- colorRampPalette(brewer.pal(5, "YlGnBu"), space="Lab")
R0Vcplot <- ggplot(combdata, aes(x = R0, y = Vc, colour=Clustering)) +
    geom_line(size=1.25) +
    scale_color_manual(values=rev(c(myPalette(5)[-1]))) +
    scale_x_continuous(trans="log", breaks=c(2, 5, 10, 15, 20), limits=c(2, 22), expand=c(0,0)) +
    scale_y_continuous(trans="probit", expand=c(0,0)) + #, breaks=seq(.5,1,.1), limits=c(.5,1)) +
    facet_grid(. ~ gx, labeller=label_parsed) +
    xlab(bquote(italic('R'[0]))) + ylab(bquote(italic("V"[c]^'*'))) +
    theme_bw() + theme(axis.text=element_text(size=9), axis.title=element_text(size=14,face="bold"))
#R0Vcplot <- R0Vcplot + coord_fixed(ratio = 1.25, xlim = NULL, ylim = NULL)
print(R0Vcplot)






#******************************************************************
#******************************************************************
####                     FIGURE 2B                             ####
#   - Phi vs Vc by R0                                             #
#******************************************************************
#******************************************************************




#---  setup Data  -------------------------------------------------

phi <- seq(1, 2, by=.005)
Vc <- seq(.8, 1, by=.005)
R0 <- seq(0, 20, by=1)

{
    vc_phi <- matrix(NA, length(Vc), length(phi))
    for (i in 1:length(Vc)){
        vc_phi[i,] <- sapply(phi, calc.R0.from.Vc, VE=.99, Vc=Vc[i])
    }
    row.names(vc_phi) <- Vc
    colnames(vc_phi)  <- phi

    vc_phi_long <- melt(vc_phi)
    colnames(vc_phi_long) <- c("Vc", "phi", "R0")
    #vc_phi_long <- vc_phi_long[which(vc_phi_long$R0<20),]
    vc_phi_long$R0[vc_phi_long$R0>=20] <- 20
    vc_phi_long <- vc_phi_long[vc_phi_long$Vc>=.80,]
    vc_phi_long$R0 <- round(vc_phi_long$R0, 2)
    vc_phi_long$phi <- round(vc_phi_long$phi, 2)
    #vc_phi_long$Vc <- round(vc_phi_long$Vc, 2)
    vc_phi_long <- unique(vc_phi_long)

    phi = sort(unique(vc_phi_long$phi))
    R0 <- sort(unique(vc_phi_long$R0))
}

#..............................................................
# Specific Disease Lines  ####

# Polio & Rubella (R0=6)
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2393518/
# http://www.who.int/immunization/polio_grad_opv_effectiveness.pdf
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2870334/pdf/16274497.pdf
vc_Polio <- sapply(phi, calc.Vc.analytic.phi, VE=1, R0=6)

# RSV (R0=3)
# http://journals.plos.org.ezp.welch.jhmi.edu/ploscompbiol/article?id=10.1371/journal.pcbi.1005133
vc_RSV <- sapply(phi, calc.Vc.analytic.phi, VE=1, R0=3)

# Cholera (R0=2.6)
# http://www.pnas.org/content/pnas/108/17/7081.full.pdf
vc_Cholera <- sapply(phi, calc.Vc.analytic.phi, VE=1, R0=2.6)


#vc_pertussis <- sapply(phi, calc.Vc.analytic.phi, VE=1, R0=10.6)
vc_pertussis <- sapply(phi, calc.Vc.analytic.phi, VE=1, R0=10.6)

# http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0083850
vc_measles <- sapply(phi, calc.Vc.analytic.phi, VE=1, R0=15)

# ://academic.oup.com/aje/article/164/10/936/162511
vc_Mumps <- sapply(phi, calc.Vc.analytic.phi, VE=1, R0=7.68)


vc_rsv <- data.frame(vc=vc_RSV, phi=phi, R0=3)
vc_chol <- data.frame(vc=vc_Cholera, phi=phi, R0=2.6)
vc_meas <- data.frame(vc=vc_measles, phi=phi, R0=15)
vc_polio <- data.frame(vc=vc_Polio, phi=phi, R0=6)
vc_mumps <- data.frame(vc=vc_Mumps, phi=phi, R0=7.68)

vc_rsv <- vc_rsv[vc_rsv$vc>=0.8,]
vc_chol <- vc_chol[vc_chol$vc>=0.8,]




#******************************************************************
# PLOT Figure 2B  ####

#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
#myPalette <- colorRampPalette(brewer.pal(9, "BuPu"), space="Lab")
#myPalette <- colorRampPalette((brewer.pal(9, "PuRd")), space="Lab")
myPalette <- colorRampPalette((brewer.pal(9, "YlOrRd")), space="Lab")

vc_htmp <- ggplot(vc_phi_long, aes(x = phi, y = Vc, z = R0))
vc_htmp <- vc_htmp + geom_tile(aes(fill = R0))
vc_htmp <- vc_htmp + scale_fill_gradientn(colours = myPalette(100),
                                          name=expression(italic(R[0]))) + theme_bw() +
    xlab(expression(italic(phi))) + ylab(expression(italic(V[c]^'*')))  +
    theme(panel.grid.minor = element_line(size=0.5), axis.text=element_text(size=10),
          axis.title=element_text(size=14,face="bold"),
          legend.key.size=unit(1,'cm'), legend.text=element_text(size=rel(1)), legend.title=element_text(size=rel(1.5))) +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
vc_htmp <- vc_htmp + geom_smooth(data=vc_meas,  aes(x=phi, y=vc), method='loess', colour="grey45", size=1)
vc_htmp <- vc_htmp + geom_smooth(data=vc_polio, aes(x=phi, y=vc), method='loess', colour="grey45", size=1)
vc_htmp <- vc_htmp + geom_smooth(data=vc_mumps, aes(x=phi, y=vc), method='loess', colour="grey45", size=1)

print(vc_htmp)



###

vc_htmp <- ggplot(vc_phi_long, aes(x = phi, y = Vc, z = R0))
vc_htmp <- vc_htmp + geom_tile(aes(fill = R0))
vc_htmp <- vc_htmp + scale_fill_gradientn(colours = myPalette(1000)[cumsum(c(1, rev(1:44)))],
                                          name=expression(italic(R[0]))) + theme_bw() +
    xlab(expression(italic(phi))) + ylab(expression(italic(V[c]^'*')))  +
    theme(panel.grid.minor = element_line(size=0.5), axis.text=element_text(size=10),
          axis.title=element_text(size=14,face="bold"),
          legend.key.size=unit(1,'cm'), legend.text=element_text(size=rel(1)), legend.title=element_text(size=rel(1.5))) +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
vc_htmp <- vc_htmp + geom_smooth(data=vc_meas,  aes(x=phi, y=vc), method='loess', colour="grey45", size=1)
vc_htmp <- vc_htmp + geom_smooth(data=vc_polio, aes(x=phi, y=vc), method='loess', colour="grey45", size=1)
vc_htmp <- vc_htmp + geom_smooth(data=vc_mumps, aes(x=phi, y=vc), method='loess', colour="grey45", size=1)

print(vc_htmp)


# For spreading colors:
cumsum(c(1, rev(1:44)))





#*********************************************************************
# Ic difference for high clustering

c(calc.Vc.analytic.phi(1, 7.7, 1), calc.Vc.analytic.phi(1, 7.7, 2))
calc.Vc.analytic.phi(1, 7.7, 2) - calc.Vc.analytic.phi(1, 7.7, 1)

c(calc.Vc.analytic.phi(1, 6, 1), calc.Vc.analytic.phi(1, 6, 2))
calc.Vc.analytic.phi(1, 6, 2) - calc.Vc.analytic.phi(1, 6, 1)

c(calc.Vc.analytic.phi(1, 2.6, 1), calc.Vc.analytic.phi(1, 2.6, 2))
calc.Vc.analytic.phi(1, 2.6, 2) - calc.Vc.analytic.phi(1, 2.6, 1)

#******************************************************************






#*********************************************************************
#*********************************************************************
#   Figure 2B - Reduced to Measles    #######
#   - R0 vs Vc for all levels of clustering where R=1
#*********************************************************************

#******************************************************************
vc_htmp2b <- ggplot(vc_phi_long, aes(x = phi, y = Vc, z = R0))
vc_htmp2b <- vc_htmp2b + geom_tile(aes(fill = R0))
vc_htmp2b <- vc_htmp2b + scale_fill_gradientn(colours = myPalette(100), name=expression(italic(R[0]))) + theme_bw() +
    xlab(expression(italic(phi))) + ylab(expression(italic(V[c]^'*')))  +
    theme(panel.grid.minor = element_line(size=0.5), axis.text=element_text(size=10), axis.title=element_text(size=14,face="bold")) +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0), limits=c(.9, 1))
# theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(), axis.ticks = element_blank(),
#       panel.background=element_rect(fill='white', colour='white'))
vc_htmp2b <- vc_htmp2b + geom_smooth(data=vc_meas, aes(x=phi, y=vc), colour="grey45", size=1)
print(vc_htmp2b)







# COMBINE INTO 1 PLOT -----------------------------------------------------

# require(gridExtra)
#

A <- R0Vcplot
B <- vc_htmp


devtools::install_github("baptiste/egg")


grid.newpage() # requires package "grid"
grid.draw(egg::ggarrange(  A + coord_fixed(ratio = 1/3),
                      B + coord_fixed(ratio = 1/.2) +
                          annotate("text", x=1.10, y=.856, label="Polio~and~Rubella~(italic(R)[0]==6.0)", angle=28, parse = TRUE, hjust = 0) +
                          annotate("text", x=1.10, y=.888, label="Mumps~(italic(R)[0]==7.7)", angle=24.5, parse = TRUE, hjust = 0) +
                          annotate("text", x=1.10, y=.944, label="Measles~(italic(R)[0]==15)", angle=14.25, parse = TRUE, hjust = 0)
                      , ncol = 1, labels=c("A","B")))















#=================================================================================================
# SUPPLEMENTAL FIGURES
#=================================================================================================

combdata <- rbind(cbind(gx=rep("A", nrow(Adata)), Adata), cbind(gx=rep("B", nrow(Bdata)), Bdata),  cbind(gx=rep("C", nrow(Cdata)), Cdata))
combdata$gx <- as.factor(combdata$gx)
levels(combdata$gx) <- list('italic(g[A](r))'='A',
                            'italic(g[B](r))'='B',
                            'italic(g[C](r))'='C')
combdata$Clustering <- factor(combdata$Clustering, levels=rev(levels(combdata$Clustering)))

R0plot_supp1 <- ggplot(combdata, aes(x = R0, y = Vc, colour=Clustering)) +
    geom_line(size=1.25) +
    scale_color_manual(values=rev(c("blue", "gold1", "darkorange1", "red3"))) +
    scale_x_continuous(breaks=c(2, 5, 10, 15, 20), limits=c(2, 22), expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    facet_grid(. ~ gx, labeller=label_parsed)
R0plot_supp1 <- R0plot_supp1 + xlab(bquote(italic('R'[0]))) + ylab(bquote(italic("V"[c]^'*')))
R0plot_supp1 <- R0plot_supp1 + theme_bw() + theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
R0plot_supp1 <- R0plot_supp1 + coord_fixed(ratio = 44, xlim = NULL, ylim = NULL)
print(R0plot_supp1)


R0plot_supp2 <- ggplot(combdata, aes(x = R0, y = Vc, colour=Clustering)) +
    geom_line(size=1.25) +
    scale_color_manual(values=rev(c("blue", "gold1", "darkorange1", "red3"))) +
    scale_x_continuous(breaks=c(2, 5, 10, 15, 20), limits=c(2, 22), expand=c(0,0)) +
    scale_y_continuous(trans="probit", expand=c(0,0)) +
    facet_grid(. ~ gx, labeller=label_parsed)
R0plot_supp2 <- R0plot_supp2 + xlab(bquote(italic('R'[0]))) + ylab(bquote(italic("V"[c]^'*')))
R0plot_supp2 <- R0plot_supp2 + theme_bw() + theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
R0plot_supp2 <- R0plot_supp2 + coord_fixed(ratio = 11, xlim = NULL, ylim = NULL)
print(R0plot_supp2)






#======================================================================================
#---  Phi vs R0 with color=Vc  ---------------------------------------------------------------------

phi <- seq(1, 5, by=.01)
R0 <- seq(0, 20, by=.01)
vc_phi <- matrix(NA, length(R0), length(phi))
for (i in 1:length(R0)){
    vc_phi[i,] <- sapply(phi, calc.Vc.analytic.phi, VE=1, R0=R0[i])
}
row.names(vc_phi) <- R0
colnames(vc_phi)  <- phi

vc_phi_long <- melt(vc_phi)
colnames(vc_phi_long) <- c("R0", "phi", "Vc")
vc_phi_long <- vc_phi_long[which(vc_phi_long$R0<25),]
head(vc_phi_long)
vc_phi_long <- vc_phi_long[which(vc_phi_long$Vc>0),]
head(vc_phi_long)
vc_phi_long <- vc_phi_long[which(vc_phi_long$R0>=5),]


vc_htmp <- ggplot(vc_phi_long, aes(x = R0, y = phi, fill = Vc)) + geom_tile()
vc_htmp <- vc_htmp + scale_fill_gradientn(colours = myPalette(100)) + theme_bw()
vc_htmp <- vc_htmp + ylab('phi') + xlab('R0')
vc_htmp <- vc_htmp + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.ticks = element_blank())
vc_htmp <- vc_htmp + ggtitle("Vaccine Thresholds: R0 vs Phi")

print(vc_htmp)

