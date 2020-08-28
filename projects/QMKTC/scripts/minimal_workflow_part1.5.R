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

load("results/cage_data_2.RData")

###################################################
###Interquantile width can also be visualized in a gene-like representation in the UCSC genome browser by exporting the data into a BED file
###################################################
exportToBed(object = myCAGEset, what = "tagClusters", qLow = 0.1, qUp = 0.9, oneFile = TRUE)
#Error in exportToBed(object = myCAGEset, what = "tagClusters", qLow = 0.1,  : 
#                       No data for given quantile positions! Run 'quantilePositions()' function for desired quantiles first, or omit 'qLow' and 'qUp' parameters to use start and end coordinates instead!

###################################################
### plot interquantile width to globally compare the promoter width across different samples
###################################################
plotInterquantileWidth(myCAGEset, clusters = "tagClusters", 
                      tpmThreshold = 3, qLow = 0.1, qUp = 0.9)
#Error in exportToBed(object = myCAGEset, what = "tagClusters", qLow = 0.1,  : 
#No data for given quantile positions! Run 'quantilePositions()' function for desired quantiles first, or omit 'qLow' and 'qUp' parameters to use start and end coordinates instead!


###################################################
### Creating consensus promoters across samples
###################################################
aggregateTagClusters(myCAGEset, tpmThreshold = 5, 
                     qLow = 0.1, qUp = 0.9, maxDist = 100)
#Error in `colnames<-`(`*tmp*`, value = "cluster") : 
#  attempt to set 'colnames' on an object with less than two dimensions
#tried this without qLow und qUp
# aggregateTagClusters(myCAGEset, tpmThreshold = 5, maxDist = 100)


###################################################
### Final set of consensus clusters can be retrieved by
###################################################
consensusCl <- consensusClusters(myCAGEset)
#which will return genomic coordinates and sum of CAGE signal across all samples for each consensus cluster (the tpm column).
head(consensusCl)
dim(consensusCl)   #

#Analogously to tag clusters, analysis of promoter width can be performed for consensus clusters as well, using the same cumulativeCTSSdistribution, quantilePositions and plotInterquantileWidth functions as described above, but by setting the clus- ters parameter to "consensusClusters".

consensusCl.tpm <- consensusClustersTpm(myCAGEset)
dim(consensusCl.tpm)
write.table(consensusCl.tpm, "results/consensus.clusters.normalized.tpm.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,  col.names = NA, qmethod = c("escape", "double"))
write.table(consensusCl,"results/consensus.clusters.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,  col.names = T, qmethod = c("escape", "double"))
rm(consensusCl,consensusCl.tpm,ctss,corr.m)

save.image("results/cage_data_2.1.RData")

system("cp -r results results_part5")
#end
