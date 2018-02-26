###########################################
### Spatial Clustering - Table 3   #####
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


