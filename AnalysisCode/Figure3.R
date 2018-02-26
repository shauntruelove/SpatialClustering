######################################################################################
### Clustering Paper 1 - Figure 3 (Heatmap-type plot of tau(r) vs 1-G(r), with Vc color)
######################################################################################


source("./R/clust_estimates_source.R")

library(SpatialClusteringPkg)

library(ggplot2); library(reshape2); library(RColorBrewer); library(grid);

library(scales)


## load Contact Parameters ####
# load(file="./data/g_params.rda")


#======================================================================================
#---  Setup Data  ---------------------------------------------------------------------

VE=1
upper.lim=1000
lambda=0.5

theta <- c(.25, .5, 1)
tau <- tau.param <- NULL

for (i in 1:length(theta)){
    tau.tmp <- calc.tau.analytic(upper.lim=upper.lim, r.delim=.1, theta=theta[i], lambda=lambda)
    tau <- rbind(tau, tau.tmp)
    tau.param.tmp <- get.tau.params(upper.lim=upper.lim, r.delim=.1, tau.level=1+theta[i], lambda=lambda)
    tau.param <- rbind(tau.param, unlist(tau.param.tmp))
}
colnames(tau) <- seq(.1, ncol(tau)*.1, .1)
row.names(tau) <- row.names(tau.param) <- theta
colnames(tau.param) <- c('theta', 'lambda', 'psi')

gamma.shape <- gamma.scale <- c(seq(.1,1,.1), seq(1.5, 5, .5), 6:10)
gamma.rate <- 1/gamma.scale

r <- seq(0, 100, .25)

vc.mat <- phi.mat <- matrix(NA, nrow=length(gamma.shape)*length(gamma.scale), ncol=length(theta))
colnames(vc.mat) <- colnames(phi.mat) <- theta
dist.list <- list(NULL)
tau.vals.list <- list(NULL)


p <- seq(0, 1, .001)
for (j in 1:length(theta)){

    count=0

    dist.mat <- tau.vals.mat <- NULL

    for (i in 1:length(gamma.shape)){
        for (h in 1:length(gamma.scale)){
            count=count+1
            dists <- qgamma(p, shape=gamma.shape[i], scale=gamma.scale[h])
            tau.vals <- calc.tau.from.pars(dists, theta=tau.param[j,'theta'], lambda=0.5, psi=tau.param[j,'psi'])

            phi.mat[count, j] <- calc.phi.cf(theta=tau.param[j,'theta'], lambda=0.5, psi=tau.param[j,'psi'], alpha=gamma.shape[i], beta=gamma.rate[h])
            vc.mat[count, j] <- 1-1/(15*phi.mat[count, j])

            dist.mat <- rbind(dist.mat, dists)
            tau.vals.mat <- rbind(tau.vals.mat, tau.vals)
        }
    }
    dist.list[[j]] <- dist.mat
    tau.vals.list[[j]] <- tau.vals.mat
}


# Curve Data for Vc Data
curve.dat <- NULL
for (j in 1:length(theta)){
    dat.tmp <- tau.vals.list[[j]]
    colnames(dat.tmp) <- p
    row.names(dat.tmp) <- 1:nrow(dat.tmp)
    dat.tmp <- melt(dat.tmp)
    colnames(dat.tmp) <- c('curve.num.clust', 'p','tau')
    dat.tmp$vc <- vc.mat[,j]
    dat.tmp$phi <- phi.mat[,j]
    dat.tmp$theta <- theta[j]
    curve.dat <- rbind(curve.dat, dat.tmp)
}
curve.dat$curve <- 1:nrow(curve.dat)
curve.dat$curveID <- paste0(curve.dat$theta,'-',curve.dat$curve.num.clust)
rm(dat.tmp)



#======================================================================================
#---  PLOT  ---------------------------------------------------------------------

library(extrafont)
#font_import()
loadfonts(device = "win")

#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
#myPalette <- colorRampPalette((brewer.pal(9, "YlGnBu")), space="Lab")
#myPalette <- colorRampPalette(rev(brewer.pal(11, "RdYlGn")), space="Lab")
myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")), space="Lab")

curve.dat <- curve.dat[rev(order(curve.dat$theta)),]
curve.dat$curveID <- factor(curve.dat$curveID, levels=(unique(curve.dat$curveID)))


# Subset for Quick Testing
curve.dat.tmp <- curve.dat[curve.dat$theta==.25,]
curve.dat.tmp <- curve.dat[curve.dat$curve.num.clust %in% sample(unique(curve.dat$curve.num.clust), 100),]
curve.dat.tmp$curveID <- factor(curve.dat.tmp$curveID, levels=(unique(curve.dat.tmp$curveID)))

plot1 <- ggplot(data=curve.dat.tmp, aes(x=p, y=tau, color = vc, group = curveID)) + geom_line(alpha=.5, size=1.5) +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
    scale_color_gradientn(name=expression(italic(I[c]^'*')), colours = myPalette(100)) + coord_fixed(ratio = 1, xlim = NULL, ylim = NULL)
plot1 + theme_bw() + theme(text=element_text(family="serif")) +
    theme(text=element_text(family="serif"), legend.position="bottom", legend.key.size=unit(1,'cm'), legend.text=element_text(size=rel(1)), legend.title=element_text(size=rel(1.5))) +
    ylab(expression(italic(tau~(r)))) + xlab(expression(1~-~italic(G(r))))



# Plot Figure 3 #########


pdf(file="./Manuscripts/Clustering Paper 1/Figures/Figure3_GrYlRd.pdf")

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdYlGn")), space="Lab")
plot1 <- ggplot(data=curve.dat, aes(x=p, y=tau, color = vc, group = curveID)) + geom_line(alpha=.25, size=1.5) +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
    scale_color_gradientn(name=expression(italic(I[c]^'*')), colours = myPalette(100)) + coord_fixed(ratio = 1, xlim = NULL, ylim = NULL)
plot1 + theme_bw() + theme(text=element_text(family="serif")) +
    theme(text=element_text(family="serif"), legend.position="bottom", legend.key.size=unit(1,'cm'), legend.text=element_text(size=rel(1)), legend.title=element_text(size=rel(1.5))) +
    ylab(expression(italic(tau~(r)))) + xlab(expression(1~-~italic(G(r))))

dev.off()


pdf(file="../../Manuscripts/Clustering Paper 1/Figures/Figure3_BuRd.pdf")

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")), space="Lab")
plot1 <- ggplot(data=curve.dat, aes(x=p, y=tau, color = vc, group = curveID)) + geom_line(alpha=.10, size=1.5) +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
    scale_color_gradientn(name=expression(italic(I[c]^'*')), colours = myPalette(100)) + coord_fixed(ratio = 1, xlim = NULL, ylim = NULL)
plot1 + theme_bw() + theme(text=element_text(family="serif")) +
    theme(text=element_text(family="serif"), legend.position="bottom", legend.key.size=unit(1,'cm'), legend.text=element_text(size=rel(1)), legend.title=element_text(size=rel(1.5))) +
    ylab(expression(italic(tau~(r)))) + xlab(expression(1~-~italic(G(r))))

dev.off()


pdf(file="./Manuscripts/Clustering Paper 1/Figures/Figure3_BuGnYl.pdf")

myPalette <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")), space="Lab")
myPalette <- colorRampPalette(c("darkblue", "turquoise","orange1"), space="Lab")

plot1 <- ggplot(data=curve.dat, aes(x=p, y=tau, color = vc, group = curveID)) + geom_line(alpha=.25, size=1.5) +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
    scale_color_gradientn(name=expression(italic(I[c]^'*')), colours = myPalette(100)) + coord_fixed(ratio = 1, xlim = NULL, ylim = NULL)
plot1 + theme_bw() + theme(text=element_text(family="serif")) +
    theme(text=element_text(family="serif"), legend.position="bottom", legend.key.size=unit(1,'cm'), legend.text=element_text(size=rel(1)), legend.title=element_text(size=rel(1.5))) +
    ylab(expression(italic(tau~(r)))) + xlab(expression(1~-~italic(G(r))))

dev.off()









#======================================================================================
# EXTRA PLOTS ####


# Phi Plot
plot1 <- ggplot(data=curve.dat, aes(x=p, y=tau, color = phi, group = phi)) + geom_line(alpha=.1, size=1.5) +
    scale_color_gradientn(name=expression(italic(phi)), colours = myPalette(100)) + coord_fixed(ratio = 1, xlim = NULL, ylim = NULL)
plot1 + theme_bw() + ylab(expression(italic(tau~(r)))) + xlab(expression(1~-~italic(G(r))))




#============================================================================
# AUC

#curve.dat$phi <- 1/((1-curve.dat$vc)*15)
curve.info <- subset(curve.dat,!duplicated(curve.dat$curveID))

curve.ids <- unique(curve.dat$curveID)
auc.vals <- numeric(length(curve.ids))
for (i in 1:length(curve.ids)){
    tmp <- curve.dat[which(curve.dat$curveID==curve.ids[i]),]
    auc.vals[i] <- sum(0.001*(tmp$tau-1))
}

auc.and.phi <- data.frame(auc=auc.vals, phi.auc=auc.vals+1, phi=curve.info$phi, theta=curve.info$theta)
cor.auc_phi <- cor(auc.and.phi$auc, auc.and.phi$phi); print(cor.auc_phi)

ggplot(data=auc.and.phi, aes(x=auc, y=phi, col=as.factor(theta))) + geom_point() +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL) +
    ylab(expression(phi)) + xlab('AUC') + guides(col=guide_legend(title=expression(theta)))

auc.and.phi$theta.lab <- factor(auc.and.phi$theta, labels=paste0('theta~"="~',unique(auc.and.phi$theta)))

ggplot(data=auc.and.phi, aes(x=auc, y=phi)) + geom_point() +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL) +
    ylab(expression(phi)) + xlab('AUC') + guides(col=guide_legend(title=expression(theta))) +
    facet_grid(. ~ theta.lab, labeller=label_parsed)



for (i in theta){
    dat.tmp <- auc.and.phi[which(auc.and.phi$theta==i),]
    print(cor(dat.tmp$auc, dat.tmp$phi))
}

