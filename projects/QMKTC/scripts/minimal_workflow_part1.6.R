#clean up
rm(list = ls(all = TRUE)) # clear all variables
graphics.off()
#set path to location where script starts, for portability
path <- getwd()
setwd(path)
path

library(CAGEr)
library(reshape)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
library(biomaRt)

load("results/cage_data_2.1.RData")


#promoter width consensus clusters
cumulativeCTSSdistribution(myCAGEset, clusters = "consensusClusters",useMulticore =T,nrCores = 4)
#error
#Error in names(clusters.cumsum) <- n : 
#  'names' attribute [37588] must be the same length as the vector [27971]
quantilePositions(myCAGEset,clusters = "consensusClusters",qLow = 0.1,qUp = 0.9,useMulticore =T,nrCores = 4) 
#error
#Getting positions of quantiles within clusters...
#-> QMKTC061AK
#Error in rep(coors$start[which(coors$tpm > 0)], ncol(cluster.q)) : 
#  invalid 'times' argument


###################################################
### plot interquantile width to globally compare the promoter width across different samples
###################################################
plotInterquantileWidth(myCAGEset, clusters = "consensusClusters", 
                       tpmThreshold = 3, qLow = 0.1, qUp = 0.9)
##Error in plotInterquantileWidth(myCAGEset, clusters = "consensusClusters",  : 
##                                  No data for given quantile positions! Run 'quantilePositions()' function for desired quantiles first!
### save data
save.image("results/cage_data_3.RData")
###

system("cp -r results results_part6")
#end

