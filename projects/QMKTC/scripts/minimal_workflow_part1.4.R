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

load("results/cage_data_1.1.RData")

###################################################
### Exporting CAGE signal to bedGraph or bigWig to visualize in UCSC Genome browser
###################################################
#exportCTSStoBedGraph(myCAGEset, values = "normalized", oneFile = TRUE,format = "bedGraph")
#exportCTSStoBedGraph(myCAGEset, values = "normalized", oneFile = FALSE,format = "bedGraph")

#save image now here:
#save.image("results/cage_data_1.2.RData")

###################################################
### Clustering
#clustering using 20bp as a maximal allowed distance between two neighbouring TSSs
#keepSingletonsAbove is useful to prevent removing highly supported singleton tag clusters
###################################################
# simple distance-based clustering
clusterCTSS(object = myCAGEset, threshold = 2, thresholdIsTpm = TRUE, 
            nrPassThreshold = 1, method = "distclu", maxDist = 20, 
            removeSingletons = TRUE, keepSingletonsAbove = 5,useMulticore = TRUE,  nrCores = 4)

myCAGEset
#Note: The clusterCTSS function creates a set of clusters for each sample separately
# export tag clusters (TCs) per sample into BED file for visualisation # in the genome browser
exportToBed(object = myCAGEset, what = "tagClusters", oneFile = TRUE)


###################################################
### Accessing promoter width
#the width of the promoter is an important characteristic that distinguishes different functional classes of promoters
#Width of every tag cluster is calculated as following:
#1. Cumulative distribution of CAGE signal along the cluster is calculated.
#2. Positions of two selected quantiles are determined
#3. Promoter width is defined as the distance (in base pairs) between the two quantiles.
###################################################
cumulativeCTSSdistribution(myCAGEset, clusters = "tagClusters" ,useMulticore = TRUE,  nrCores = 4)
##gives error, Error in names(clusters.cumsum) <- n : 
#'names' attribute [83044] must be the same length as the vector [50664]
quantilePositions(myCAGEset, clusters = "tagClusters", qLow = 0.1, qUp = 0.9 ,useMulticore = TRUE,  nrCores = 4)
##gives error

save.image("results/cage_data_2.RData")

system("cp -r results results_part4")
#end
