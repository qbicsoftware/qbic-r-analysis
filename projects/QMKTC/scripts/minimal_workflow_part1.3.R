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

load("results/cage_data_1.RData")

###################################################
### Genomic coordinates of all TSSs and numbers of supporting CAGE tags in every input sample can be retrieved by typing
###################################################
library(genefilter)
f1 <- kOverA(6, 1)
ffun <- filterfun(f1)
wh1 <- genefilter(myCAGEset@tagCountMatrix, ffun)
myCAGEset@tagCountMatrix <- myCAGEset@tagCountMatrix[wh1,]
myCAGEset@CTSScoordinates <- myCAGEset@CTSScoordinates[wh1,]

dim(myCAGEset@CTSScoordinates)
dim(myCAGEset@tagCountMatrix)

ctss <- CTSStagCount(myCAGEset)
class(ctss);head(ctss);dim(ctss)


###################################################
### check sample labels and colour assigned
###################################################
sampleLabels(myCAGEset)
myCAGEset@inputFiles


# ###################################################
# ### correlation between the samples
# ###################################################
corr.m <- cor(ctss[,c(4:ncol(ctss))], method = "pearson")
write.table(corr.m,"results/correlation.matrix.txt",sep=" \t")


###################################################
### check library sizes
## Library sizes (number of total sequenced tags) of individual experiments differ, thus normalization is required to make them comparable
###################################################
librarySizes(myCAGEset)
libs <- data.frame(libsizes=librarySizes(myCAGEset),QBiC.Code=sampleLabels(myCAGEset))
#libs$ID <- row.names(libs)
#names(libs) <- c("libsizes","QBiC.Code")
libs <- merge(libs,meta,by="QBiC.Code")
#plot
pdf("results/libsizes.pdf")
par(mar=c(3,3,3,3))
ggplot(data=libs, aes(x=QBiC.Code, y=libsizes,fill=Condition..treatment)) + geom_bar(stat="identity",position=position_dodge(),color="black") +
  #scale_fill_manual(values=c("#E69F00", "#999999")) +
  facet_grid(~Condition..cell_line+Condition..time+Condition..genotype,scales="free",space = "free") +
  ylab('Library sizes (number of total sequenced tags)') +
  xlab('Samples') +
  ggtitle('') +
  theme(text = element_text(size=6),
        axis.text.x = element_text(angle=60, vjust=1,hjust=1)) #vjust and hjust adjust the vertical and horizontal justification of the labels, which is often useful
dev.off()
write.table(libs, "results/librarySizes.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,  col.names = NA, qmethod = c("escape", "double"))



###################################################
#Normalization
### does dataset follow power- law distribution?
#It has been shown that many CAGE datasets follow a power- law distribution (Balwierz et al. (2009))
#To check whether our CAGE datasets follow power-law distribution and in which range of values, we can use the plotReverseCumulatives function:
###################################################
plotReverseCumulatives(myCAGEset, values="raw", fitInRange = c(5, 10000), onePlot = TRUE)

###################################################
### normalization using powerLaw here
###################################################
#use library(multicore) here in the future??
normalizeTagCount(myCAGEset, method = "none", 
                  fitInRange = c(5, 10000), alpha = 1.17, T = 1e+07)
# method "none"would allow to carry on the raw tag counts

#plot reverse cumulatives with normalized tags
plotReverseCumulatives(myCAGEset, values="normalized", fitInRange = c(5, 10000), onePlot = TRUE)
save.image("results/cage_data_1.1.RData")

system("cp -r results results_part3")
#end
