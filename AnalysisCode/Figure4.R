######################################################################################
### Clustering Paper 1 - Figure 3 (Heatmap-type plot of tau(r) vs 1-G(r), with Vc color)
######################################################################################

#source('./R/multiplot.R')
#source("./R/clust_estimates_source.R")

## load Contact Parameters ####
# load(file="./data/g_params.rda")

#install.packages("SpatialClusteringPkg_0.1.0.tar.gz", repos = NULL, type="source")
library(SpatialClusteringPkg)
library(ggplot2)
library(RColorBrewer)
library(reshape2)


######################################################################################
### Clustering Paper 1 - Figure 4 (Heatmap-type plot of phi vs vacc, with RR outbreak color)
######################################################################################

recalc.data <- FALSE




#======================================================================================
#---  Setup Data  ---------------------------------------------------------------------

if (recalc.data) {
  
  VE=1
  phi <- seq(1, 2, .025)
  #vacc <- c(seq(.75,.85,.05), seq(.81,1,.01))
  vacc <- seq(.75,1,.01)
  #vacc <- vacc[-length(vacc)]
  
  
  
  # Eff. Reprod Number for measles
  R0 <- 15
  R15 <- R0 * (t(t((1-vacc*VE))) %*% t(phi))
  colnames(R15) <- phi
  row.names(R15) <- vacc
  
  
  # Outbreak Probability 
  source("./Source/StochasticSIR_source.R")
  
  prob.mat <- matrix(NA, nrow=length(vacc), ncol=length(phi))
  colnames(prob.mat) <- phi
  row.names(prob.mat) <- vacc
  for (i in 1:length(vacc)){
    for (j in 1:length(phi)){
      prob.mat[i,j] <- sum((simulateFinalSizes(n=10000,evt.driven=FALSE, beta=R15[i,j], gamma=1, 
                                               initial.state = c(S = 999, I = 1, R = 0), freq.dependent = TRUE, step.size = .1)/1000)>=.05)/10000
    }
  }
  
  prob.unclust <- prob.mat[,1]
  RR.mat <- prob.mat / prob.unclust
  
  
  
  R15.dat <- list(R0=R0, VE=VE, phi=phi, vacc=vacc, R=R15, prob.mat=prob.mat, prob.unclust=prob.unclust, RR.mat=RR.mat)
  
  # Save for future use --> Takes a long time to calculate:
  save(R15.dat, file='./Manuscripts/Clustering Paper 1/Figures/Figure4_R15_data.RData')
  
  load(file='./Manuscripts/Clustering Paper 1/Figures/Figure4_R15_data.RData') # Loads R15.dat
  
  
  
  prob.vacc80 <- prob.mat[1,]
  RR.mat.vacc80 <- t(t(prob.mat) / prob.vacc80)
  
  RR.vacc.mat <- matrix(NA, nrow=nrow(RR.mat), ncol=ncol(RR.mat))
  for (i in 1:(nrow(prob.mat)-1)){
    prob.vacc0 <- prob.mat[i,]
    RR.vacc.mat[i,] <- prob.mat[i+1,] / prob.vacc0
  }
  colnames(RR.vacc.mat) <- phi
  row.names(RR.vacc.mat) <- vacc
  RR.vacc.mat <- 1/RR.vacc.mat
  
  
  
  
  #..............................................................................................
  # R0=6 #########
  
  
  # Eff. Reprod Number for measles
  R0 <- 6
  
  R6 <- R0 * (t(t((1-vacc*VE))) %*% t(phi))
  colnames(R6) <- phi
  row.names(R6) <- vacc
  
  
  # Outbreak Probability 
  source("./Source/StochasticSIR_source.R")
  
  prob.mat <- matrix(NA, nrow=length(vacc), ncol=length(phi))
  colnames(prob.mat) <- phi
  row.names(prob.mat) <- vacc
  for (i in 1:length(vacc)){
    for (j in 1:length(phi)){
      prob.mat[i,j] <- sum((simulateFinalSizes(n=10000,evt.driven=FALSE, beta=R6[i,j], gamma=1, 
                                               initial.state = c(S = 999, I = 1, R = 0), freq.dependent = TRUE, step.size = .1)/1000)>=.05)/10000
    }
  }
  
  prob.unclust <- prob.mat[,1]
  RR.mat <- prob.mat / prob.unclust
  
  R6.dat <- list(R0=R0, VE=VE, phi=phi, vacc=vacc, R=R6, prob.mat=prob.mat, prob.unclust=prob.unclust, RR.mat=RR.mat)
  
  # Save for future use --> Takes a long time to calculate:
  save(R6.dat, file='./Manuscripts/Clustering Paper 1/Figures/Figure4_R6_data.RData')
  
  load(file='./Manuscripts/Clustering Paper 1/Figures/Figure4_R6_data.RData') # loads R6.dat
  list2env(R6.dat, .GlobalEnv)
  
}








#======================================================================================
#---  R=15  PLOT  ---------------------------------------------------------------------

# Unlist Saved objects to Environment
load(file='./Manuscripts/Clustering Paper 1/Figures/Figure4_R15_data.RData') # loads R15.dat
list2env(R15.dat, .GlobalEnv)



# Probability Plot --------------------------------------------------------

prob.long <- melt(prob.mat[which(row.names(prob.mat)<=1 & row.names(prob.mat)>=.75),])
prob.long <- melt(prob.mat)

colnames(prob.long) <- c("Vacc", "Phi", "Prob")
vacc.lims <- c(.75, 1.0)
prob.lims <- c(min(prob.long$Prob), max(prob.long$Prob))

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
myPalette <- colorRampPalette((brewer.pal(9, "PuRd")), space="Lab")
myPalette <- colorRampPalette(c("lightgreen", "navy"), space="Lab")

htmp1 <- ggplot(prob.long, aes(x = Vacc, y = Phi, fill = Prob)) + geom_tile()
htmp1 <- htmp1 + scale_fill_gradientn(colours = myPalette(1000), limits=prob.lims) + 
  ylab(expression(italic(phi))) + xlab(expression(italic(v))) + 
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.ticks = element_blank()) + theme_bw() + 
  coord_fixed(ratio = ((vacc.lims[2]-vacc.lims[1])/1), xlim = vacc.lims, ylim = NULL) + ggtitle('Outbreak Probability')
print(htmp1)



# RR Plot -----------------------------------------------------------------

# Fix NAs
prob.unclust.tmp <- prob.unclust
prob.unclust.tmp[prob.unclust.tmp==0] <- min(prob.mat[prob.mat>0])/2
RR.mat <- prob.mat / prob.unclust.tmp

#RR.long <- melt(RR.mat[which(row.names(RR.mat)<=.96 & row.names(RR.mat)>=.8),])
RR.long <- melt(RR.mat)

colnames(RR.long) <- c("Vacc", "Phi", "RR")
RR.long$RR[which(RR.long$RR>=10)] <- 10
RR.long$RR[which(RR.long$RR<=1)] <- 1
RR.long$RR[which(RR.long$Vacc>.97)] <- NA

vacc.lims <- c(.75, 1)

myPalette <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")), space="Lab")
myPalette <- colorRampPalette((brewer.pal(9, "YlOrRd")), space="Lab")
myPalette <- colorRampPalette(rev(heat.colors(1000, alpha = 1)), space="Lab")
myPalette <- colorRampPalette(c("cyan3", "darkblue"), space="Lab")
myPalette <- colorRampPalette(c("khaki1", "red4"), space="Lab")

htmp2 <- ggplot(RR.long, aes(x = Vacc, y = Phi, fill = RR)) + geom_tile()
htmp2 <- htmp2 + scale_fill_gradientn(colours = myPalette(1000), breaks=c(1,3,5,7.5,10), labels=c(1,3,5,7.5,"\u2265 10"), limits=c(1,10), na.value='white') + 
  ylab(expression(italic(phi))) + xlab(expression(italic(v))) + 
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.ticks = element_blank()) + theme_bw() + 
  coord_fixed(ratio = ((vacc.lims[2]-vacc.lims[1])/1), xlim = vacc.lims, ylim = NULL) + ggtitle('Relative Outbreak Risk')
print(htmp2)



multiplot(htmp1, htmp2, cols=2)

# Create new ggplot objects so we can print all together
htmp1.15 <- htmp1
htmp2.15 <- htmp2


# Better way to print together --------------------------------------------

#devtools::install_github("baptiste/egg")

library(egg); library(grid)

pdf(file="./Manuscripts/Clustering Paper 1/Figures/Figure4_R15.pdf") 

grid.newpage()

grid.draw(ggarrange(htmp1, htmp2, ncol = 2, labels=c("A","B")))

dev.off()






# #------------------------------------
# # Log Scale
# 
# RR.mat.log <- log(RR.mat)
# 
# RR.long <- melt(RR.mat.log)
# colnames(RR.long) <- c("Vacc", "Phi", "RR")
# htmp1 <- ggplot(RR.long, aes(x = Vacc, y = Phi, fill = RR)) + geom_tile()
# htmp1 <- htmp1 + scale_fill_gradientn(colours = myPalette(100)) + theme_bw() + 
#   ylab(expression(phi)) + xlab(expression(italic(v))) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
#                        panel.grid.minor = element_blank(), axis.ticks = element_blank()) + 
#   coord_fixed(ratio = 1/2, xlim = NULL, ylim = NULL)
# print(htmp1)
# 
# 








#======================================================================================
#---  R=6  PLOT  ---------------------------------------------------------------------

# Unlist Saved objects to Environment
load(file='./Manuscripts/Clustering Paper 1/Figures/Figure4_R6_data.RData') # loads R15.dat
list2env(R6.dat, .GlobalEnv)



# Probability Plot --------------------------------------------------------

prob.long <- melt(prob.mat[which(row.names(prob.mat)<=1 & row.names(prob.mat)>=.75),])
prob.long <- melt(prob.mat)
colnames(prob.long) <- c("Vacc", "Phi", "Prob")
vacc.lims <- c(.75, 1)


myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
myPalette <- colorRampPalette(rev(brewer.pal(9, "PuRd")), space="Lab")
myPalette <- colorRampPalette((brewer.pal(9, "YlOrRd")), space="Lab")
myPalette <- colorRampPalette(c("lightgreen", "navy"), space="Lab")


htmp1 <- ggplot(prob.long, aes(x = Vacc, y = Phi, fill = Prob)) + geom_tile()
htmp1 <- htmp1 + scale_fill_gradientn(colours = myPalette(1000), limits=prob.lims) + 
  ylab(expression(italic(phi))) + xlab(expression(italic(v))) + 
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.ticks = element_blank()) + theme_bw() + 
  coord_fixed(ratio = ((vacc.lims[2]-vacc.lims[1])/1), xlim = vacc.lims, ylim = NULL) + ggtitle('Outbreak Probability')
print(htmp1)



# RR Plot ----------------------------------------------------------------------

# Fix NAs
prob.unclust.tmp <- prob.unclust
prob.unclust.tmp[prob.unclust.tmp==0] <- min(prob.mat[prob.mat>0])/2
RR.mat <- prob.mat / prob.unclust.tmp

RR.long <- melt(RR.mat[which(row.names(RR.mat)<=.96 & row.names(RR.mat)>=.8),])
RR.long <- melt(RR.mat)

colnames(RR.long) <- c("Vacc", "Phi", "RR")
RR.long$RR[which(RR.long$RR>=10)] <- 10
RR.long$RR[which(RR.long$RR<=1)] <- 1
RR.long$RR[which(RR.long$Vacc>=.94)] <- NA

# Color Palette
myPalette <- colorRampPalette((brewer.pal(9, "YlGnBu")), space="Lab")
myPalette <- colorRampPalette((brewer.pal(9, "YlOrRd")), space="Lab")
myPalette <- colorRampPalette(rev(heat.colors(1000, alpha = 1)), space="Lab")
myPalette <- colorRampPalette(c("khaki1", "red4"), space="Lab")

htmp2 <- ggplot(RR.long, aes(x = Vacc, y = Phi, fill = RR)) + geom_tile()
htmp2 <- htmp2 + scale_fill_gradientn(colours = myPalette(1000), breaks=c(1,3,5,7.5,10), labels=c(1,3,5,7.5,"\u226510"), limits=c(1,10), na.value='white') + 
  ylab(expression(italic(phi))) + xlab(expression(italic(v))) + 
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  #scale_colour_manual(values=NA) + guides(colour=guide_legend("Prob=0", override.aes=list(colour="white"))) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.ticks = element_blank()) + theme_bw() + 
  coord_fixed(ratio = ((vacc.lims[2]-vacc.lims[1])/1), xlim = vacc.lims, ylim = NULL) + ggtitle('Relative Outbreak Risk')


print(htmp2)  
multiplot(htmp1, htmp2, cols=2)


# Better way to print together --------------------------------------------

#devtools::install_github("baptiste/egg")

library(egg); library(grid)

pdf(file="./Manuscripts/Clustering Paper 1/Figures/Figure4_R6.pdf") 

grid.newpage()
grid.draw(ggarrange(htmp1, htmp2, ncol = 2, labels=c("A","B")))

dev.off()




#.......................................................................................
## LOG SCALE RR ####

# Fix NAs
prob.unclust.tmp <- prob.unclust
prob.unclust.tmp[prob.unclust.tmp==0] <- min(prob.mat[prob.mat>0])/2
RR.mat <- prob.mat / prob.unclust.tmp

#RR.long <- melt(RR.mat[which(row.names(RR.mat)<=.96 & row.names(RR.mat)>=.8),])
RR.long2 <- melt(RR.mat)
colnames(RR.long2) <- c("Vacc", "Phi", "RR")

# RR.long$RR[which(RR.long$RR>=10)] <- 10
# RR.long$RR[which(RR.long$RR<=1)] <- 1
# RR.long$RR[which(RR.long$Vacc>=.97)] <- NA

RR.long.log <- RR.long2
RR.long.log$RR <- log(RR.long2$RR)

RR.long.log$RR[which(RR.long2$RR==0)] <- NA

summary(RR.long2$RR); summary(RR.long.log$RR)


htmp2log <- ggplot(RR.long.log, aes(x = Vacc, y = Phi, fill = RR)) + geom_tile()
htmp2log <- htmp2log + scale_fill_gradientn(colours = myPalette(1000), breaks=log(c(1,3,5,7.5,10)), labels=c(1,3,5,7.5,"\u226510"), limits=log(c(1,10)), na.value='white') + 
  ylab(expression(italic(phi))) + xlab(expression(italic(v))) + 
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.ticks = element_blank()) + theme_bw() + 
  coord_fixed(ratio = ((vacc.lims[2]-vacc.lims[1])/1), xlim = vacc.lims, ylim = NULL) + ggtitle('Relative Outbreak Risk')
print(htmp2log)

#.......................................................................................












# PLOT ALL 4 TOGETHER -----------------------------------------------------


#devtools::install_github("baptiste/egg")

library(egg); library(grid)

{ pdf(file="./Manuscripts/Clustering Paper 1/Figures/Figure4_all_compact.pdf", onefile=FALSE) 
  
  grid.newpage()
  grid.draw(ggarrange(
    htmp1 + xlab(NULL) + theme(legend.position="none") + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank()),
    htmp2 + xlab(NULL) + theme(legend.position="none") + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
            axis.title.y=element_blank(), axis.text.y=element_blank()), 
    htmp1.15 + ggtitle(NULL) + 
      theme(legend.position="bottom", legend.key.size=unit(.75,'cm'), legend.text=element_text(size=rel(1)), legend.title=element_text(size=rel(1.5))),
    htmp2.15 + ggtitle(NULL) + 
      theme(legend.position="bottom", legend.key.size=unit(.75,'cm'), legend.text=element_text(size=rel(1)), legend.title=element_text(size=rel(1.5)), 
            axis.title.y=element_blank(), axis.text.y=element_blank()),
    ncol = 2))
  
  dev.off() }



# Compact Plot Flipped ----------------------------------------------------

{ pdf(file="./Manuscripts/Clustering Paper 1/Figures/Figure4_all_compact2.pdf", onefile=FALSE) 
  
  grid.newpage()
  grid.draw(ggarrange(
    
    htmp1 + xlab(NULL) + theme(legend.position="none") + ggtitle(expression(italic(R[0])==6)) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), plot.title = element_text(hjust = 0.5)),
    
    htmp1.15 + xlab(NULL) + theme(legend.position="none") + ggtitle(expression(italic(R[0])==15)) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
            axis.title.y=element_blank(), axis.text.y=element_blank(), 
            plot.title = element_text(hjust = 0.5), 
            legend.position="right", legend.key.size=unit(.75,'cm'), legend.text=element_text(size=rel(.75)), 
            legend.title=element_text(size=rel(1))),
    
    htmp2 + ggtitle(NULL) + theme(legend.position="none") +
      annotate("text", x=mean(c(.935,1)), y=1.5, colour='grey70', angle=90, label="italic(Pr(outbreak))==0", parse=TRUE) +
      theme(axis.text.x=element_text(angle=45)),
    
    htmp2.15 + ggtitle(NULL) + 
      theme(legend.position="right", legend.key.size=unit(.75,'cm'), legend.text=element_text(size=rel(.75)), legend.title=element_text(size=rel(1)), 
            axis.title.y=element_blank(), axis.text.y=element_blank(), axis.text.x=element_text(angle=45)) +
      annotate("text", x=mean(c(.975,1)), y=1.5, colour='grey70', angle=90, label="italic(Pr(outbreak))==0", parse=TRUE),
    
    ncol = 2))
  
  dev.off() }



# Compact Plot Flipped ----------------------------------------------------

{ png(file="./Manuscripts/Clustering Paper 1/Figures/Figure4_all_compact2.png", 
      width=15, height=15, units='in', res=600, type='cairo', bg='transparent') 
  
  grid.newpage()
  grid.draw(ggarrange(
    
    htmp1 + xlab(NULL) + theme(legend.position="none") + ggtitle(expression(italic(R[0])==6)) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), plot.title = element_text(hjust = 0.5)),
    
    htmp1.15 + xlab(NULL) + theme(legend.position="none") + ggtitle(expression(italic(R[0])==15)) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
            axis.title.y=element_blank(), axis.text.y=element_blank(), 
            plot.title = element_text(hjust = 0.5), 
            legend.position="right", legend.key.size=unit(.75,'cm'), legend.text=element_text(size=rel(.75)), 
            legend.title=element_text(size=rel(1))),
    
    htmp2 + ggtitle(NULL) + theme(legend.position="none") +
      annotate("text", x=mean(c(.935,1)), y=1.5, colour='grey70', angle=90, label="italic(Pr(outbreak))==0", parse=TRUE) +
      theme(axis.text.x=element_text(angle=45)),
    
    htmp2.15 + ggtitle(NULL) + 
      theme(legend.position="right", legend.key.size=unit(.75,'cm'), legend.text=element_text(size=rel(.75)), legend.title=element_text(size=rel(1)), 
            axis.title.y=element_blank(), axis.text.y=element_blank(), axis.text.x=element_text(angle=45)) +
      annotate("text", x=mean(c(.975,1)), y=1.5, colour='grey70', angle=90, label="italic(Pr(outbreak))==0", parse=TRUE),
    
    ncol = 2))
  
  dev.off() }


# Compact Plot with Boxes -------------------------------------------------

# geom_rect(aes(xmin=.915, xmax=.975, ymin=.99, ymax=2.01), fill=NA, colour='red', size=1.25, inherit.aes = FALSE)
# 
# htmp2 <- htmp2 + 
#   geom_rect(aes(xmin=.915, xmax=.975, ymin=.99, ymax=2.01), fill=NA, colour='red', size=1.25, inherit.aes = FALSE) +
#   annotate("text", x=mean(c(.975,1)), y=1.5, colour='grey65', angle=90, label="italic(Pr(Outbreak))==0", parse=TRUE)





# grid.newpage()
# grid.draw(ggarrange(
#   htmp1.15 + xlab(NULL) + theme(legend.position="none") + 
#     theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()),
#   htmp2.15 + xlab(NULL) + theme(legend.position="none") + 
#     theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
#           axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()),
#   htmp1 + ggtitle(NULL) + 
#     theme(legend.position="bottom", legend.key.size=unit(.75,'cm'), legend.text=element_text(size=rel(1)), legend.title=element_text(size=rel(1.5))),
#   htmp2 + ggtitle(NULL) + 
#     theme(legend.position="bottom", legend.key.size=unit(.75,'cm'), legend.text=element_text(size=rel(1)), legend.title=element_text(size=rel(1.5)), 
#           axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()),
#   ncol = 2))



{ pdf(file="./Manuscripts/Clustering Paper 1/Figures/Figure4_all.pdf", onefile=FALSE) 
  
  #grid.newpage()
  grid.draw(ggarrange(
    htmp1.15, htmp2.15, htmp1, htmp2,
    ncol = 2))
  
  dev.off() }




