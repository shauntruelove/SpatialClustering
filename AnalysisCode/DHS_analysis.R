######################################################################################
### Clustering Paper 1 - DHS ANALYSIS (Table 3 and Figure 5)
######################################################################################

# source("./FinalSize/Source/finalsize_source.R")
# source("./FinalSize/Source/finalsize_utilities.R")
# source("./R/sim_source.R")
#source("./R/clust_estimates_source.R")
# source("./DHS_Estimation/Source/DHSClustEst_source.R")
# source("./DHS_Estimation/Source/DHSClustEst_map_source.R")

library(SpatialClusteringPkg)
source("DHS_source.R")
source("Paper1_github/Source/DHS_source.R")



#file.edit("./DHS_Estimation/Source/DHSClustEst_source.R")
# 
# library(ggplot2)
# library(reshape2) # melt
# library(RColorBrewer) # brewer.pal
# library(rgdal)  # readOGR
# library(SpatialEpi) # latlong2grid
# library(foreign) # read.dta
# #library(maptools) # readShapeSpatial
# 

############################################################################################
###
###       FIGURE 5 - DHS Map                                 #######
###
############################################################################################

# TANZANIA DHS DATA -------------------------------------------------------

country <- 'Tanzania'
year <- 2010
country.dhs <- "TZKR63FL.DTA"
country.gps <- "TZGE61FL.dbf"
directory <- "C:/Users/Shaun/JHSPH/IDD SVN Repository/Measles/data"


#  MAP TANZANIA ------------------------------------------------------
# Import DHS Shapefile for the country. This is available from DHS at https://spatialdata.dhsprogram.com/boundaries/#view=table&countryId=AF.
#  These shapefiles provide by DHS should be used because they have the correct boundaries and region/districts for the DHS data.
tanz2010_map <- map.country.clust.gg(country=country, year=year, age.range=c(13,60),
                                        filename=country.dhs, filename.gps=country.gps, treat_NAs="unvacc",
                                        directory=directory, shapefile.dir=file.path(directory, 'DHS_boundaries'),
                                        pt.size=2, pt.transpar=.60, scale.dist=50, scale.loc='bottomright')

print(tanz2010_map)










############################################################################################
###
###       DHS Spatial Clustering - Table 3                                  #######
###
############################################################################################

library(SpatialClusteringPkg)
source("DHS_source.R")
source("Paper1_github/Source/DHS_source.R")
library(reshape2)

# source("./R/outbreak_prop_functs.R")
# source("./FinalSize/Source/finalsize_utilities.R")

write.excel <- function(x,row.names=TRUE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}


# TANZANIA DHS DATA -------------------------------------------------------

country <- 'Tanzania'
year <- 2010
country.dhs <- "TZKR63FL.DTA"
country.gps <- "TZGE61FL.dbf"
directory <- "C:/Users/Shaun/JHSPH/IDD SVN Repository/Measles/data"
iso.tmp <- maps::iso.alpha(country, 3)

### DHS DATA

lower_age <- 12
upper_age <- 60
lambda=2
boot.iter=1000
theta=.5
treat_NAs <- "remove"


r.max=c(.1,.5,1,2,4,6,9,12,16,20,25,50,100,200)
r.min <- c(0, r.max[-length(r.max)])
#r.min <- c(0, r.max[-length(r.max)])
r.mid <- (r.max+r.min)/2
n_r <- length(r.max); n_r

plot.maps = FALSE
district.taus = FALSE

# Import, Clean, and Merge DHS data using functions from "DHS_source.R"

dat <- dhs.data.import(country, data.year=year, filename=country.dhs, filename.gps=country.gps,
                       directory=directory, 
                       youngest.age=0, oldest.age=60, 
                       model.type="all", merge.gps=TRUE)
dat1 <- dat$data
# Remove points without coordinates
dat1 <- dat1[!(round(dat1$x.coord.km,5)==0 & round(dat1$y.coord.km,5)==0),]

var.labels <- dat$var.labels
dat.clust <- aggr.dhs.clusters(dat1, vaccbreaks=seq(0,1,.1),  treat_NAs=treat_NAs)
dat.clust <- dat.clust[!(dat.clust$x==0 & dat.clust$y==0),]
n.clusts <- nrow(dat.clust)
dat <- NULL


#///////////////////////////////////////////////////////////////////////////////////////////////////
# Cluster Distances & Sizes

clust.dists <- dist(dat.clust[,c('x','y')], method='euclidean')
clust.dists.mat <- as.matrix(clust.dists)
clust.dists.mat[which(clust.dists.mat==0)] <- NA
dists.vect <- as.vector(clust.dists)
max.dist <- ceiling(max(dists.vect)/100)*100






#/////////////////////////////////////////////////////////////////////////////////////////////////////////////
#                       Calculate Taus                        #########
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////


# # Original Cluster points
# colnames(dat.clust)
# clust.pop.dhs <- data.frame(as.data.frame(dat.clust[,c('vaccrate','x','y','num.indiv', 'sampwt','urban')]))
# colnames(clust.pop.dhs) <- c('vaccrate', 'x', 'y', 'num.indiv', 'sampwt','urban')
# 


#/////////////////////////////////////////////////////////////////////////////////////////////////////////////
# Clean up data sources for DHS Tau Estimate

# Aggregate individual data to cluster-based data
dat.clust <- aggr.dhs.clusters(dat1, vaccbreaks=seq(0,1,.1), treat_NAs=treat_NAs)

# Create separate datasets for each specific tau type
dat.clust.wt <- dat.clust

# make sample weight 1 for non-weighted taus
dat.clust[,'sampwt'] <- 1

dat.clust    <- DHS.to.spatialtau.clust(data=dat.clust)
dat.clust.wt <- DHS.to.spatialtau.clust(data=dat.clust.wt)


#===========================================================================================================================
# Calculate Taus - jitter, boot, etc.
# r.max=c(.1,1,2,4,6,9,12,16,20,25)
# lambda=2
# boot.iter=20
plot.dir=NA
#plot.dir=paste0('./DHS_Estimation/Results/GAM_ProbSurfaces/', DHS_code, '_TauFits')
# theta=.5


tau.data=dat.clust; country=country; data.year=year; remove.hh=FALSE; sampwt=FALSE; 
remove.outliers=remove.outliers; outlier.cutoff=outlier.cutoff;
r.max=r.max; fit.lambda=T; window.tau=FALSE; plot.fits=FALSE;
theta=theta; lambda=lambda; district=NA; district.excl=NA; plot.taus=FALSE; boot.iter=boot.iter;
plot.dir=plot.dir; save.plots=FALSE; remove_self=FALSE; blend.first.r=T


taufit       <- calc.tau.dhs.data(tau.data=dat.clust, country=country, data.year=year, remove.hh=FALSE, sampwt=FALSE, 
                                  r.max=r.max, fit.lambda=T, plot.fits=FALSE,
                                  theta=theta, lambda=lambda, plot.taus=FALSE, boot.iter=boot.iter,
                                  plot.dir=plot.dir, save.plots=FALSE, blend.first.r=T)

taufit.wt    <- calc.tau.dhs.data(tau.data=dat.clust.wt, country=country, data.year=year, remove.hh=FALSE, sampwt=TRUE, 
                                  remove.outliers=remove.outliers, outlier.cutoff=outlier.cutoff,
                                  r.max=r.max, fit.lambda=FALSE, window.tau=FALSE, plot.fits=FALSE,
                                  theta=theta, lambda=lambda, district=NA, district.excl=NA, plot.taus=FALSE, boot.iter=boot.iter,
                                  plot.dir=plot.dir, save.plots=FALSE, blend.first.r=T)

save(taufit.hh, taufit.wt.hh, r.max, file=paste0('./DHS_Estimation/',results.dir,DHS_code,'_taufits_cluster_all.RData'))



pdf(file=paste0('./DHS_Estimation/',results.dir, DHS_code, '_ParamDistribs.pdf'))

par(mfrow=c(4,2))
# Plot Fitted Parameters
hist(taufit$fits[,1], breaks=100, xlab='theta', main='Fitted Thetas')
hist(taufit$fits[,2], breaks=100, xlab='lambda', main='Fitted Lambdas')
hist(taufit.hh$fits[,1], breaks=100, xlab='theta', main='Fitted Thetas')
hist(taufit.hh$fits[,2], breaks=100, xlab='lambda', main='Fitted Lambdas')
hist(taufit.wt$fits[,1], breaks=100, xlab='theta', main='Fitted Thetas')
hist(taufit.wt$fits[,2], breaks=100, xlab='lambda', main='Fitted Lambdas')
hist(taufit.wt.hh$fits[,1], breaks=100, xlab='theta', main='Fitted Thetas')
hist(taufit.wt.hh$fits[,2], breaks=100, xlab='lambda', main='Fitted Lambdas')

dev.off()  






