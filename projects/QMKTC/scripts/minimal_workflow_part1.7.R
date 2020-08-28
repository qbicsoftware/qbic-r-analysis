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
library(rtracklayer)
library(ChIPpeakAnno)

tpm = read.table("results/consensus.clusters.normalized.tpm.txt",header = T,sep = "\t",
                 na.strings =c("","NaN"),quote=NULL,stringsAsFactors=F,dec=".",fill=TRUE,row.names = 1)
head(tpm)

cor = cor(tpm,method = "pearson")
cor = melt(cor)

pdf("correlation_consensus_cluster_tpm_data.pdf")
ggplot(data = cor, aes(x=X1, y=X2, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  ylab("") +
  xlab("") +
  theme(text = element_text(size=6),
        axis.text.x = element_text(angle=90, vjust=1,hjust=1))
dev.off()
#put matrix to file
write.table(cor, "results/correlation.matrix_consensus_cluster_tpm_data.txt",sep="\t")


#get chromosome count distribution
rm(list = ls(all = TRUE)); graphics.off()
path <- getwd(); setwd(path)
path
d1 = read.table("results/consensus.clusters.txt",header = T,sep = "\t",na.strings =c("","NaN"),quote=NULL,stringsAsFactors=F,dec=".",fill=TRUE)
d1$ID = paste("cl_",d1$consensus.cluster,sep="")
d1$consensus.cluster <- NULL
d1 = d1[,c("chr","ID")]
d2 = read.table("results/consensus.clusters.normalized.tpm.txt",header = T,sep = "\t",na.strings =c("","NaN"),quote=NULL,stringsAsFactors=F,dec=".",fill=TRUE)
d2$ID = paste("cl_",d2$X,sep="")
d2$X <- NULL
me = merge(d1,d2,by="ID")
me = melt(me,id.vars = c("ID","chr"))
me$QBiC.Code =me$variable
m <- read.table("results/metadata.txt",header = T,sep = "\t",na.strings =c("","NaN"),quote=NULL,stringsAsFactors=F,dec=".",fill=TRUE)
me <- merge(me,m,by="QBiC.Code")
head(me)

pdf("chromosome_consensus_cluster_tpm_distribution.pdf")
ggplot(data=me, aes(x=QBiC.Code, y=value,fill=Condition..treatment)) + geom_bar(stat = "identity") +
  ylab('tpm') +
  xlab('') +
  ggtitle('') +
  facet_grid(~Condition..cell_line+Condition..time+Condition..genotype,scales="free",space = "free") +
  theme(text = element_text(size=6),
        axis.text.x = element_text(angle=45, vjust=1,hjust=1)) #vjust and hjust adjust the vertical and 
dev.off()

#get plot about distribution across promoter region,introns,etc see cager vignette to plotannot()
#annotate first
rm(list = ls(all = TRUE)); graphics.off()
path <- getwd(); setwd(path)
path

setwd("Ensembl38.97/")
GTFfile = "Homo_sapiens.GRCh38.97.gtf"
GTFfile
GTF <- import.gff(GTFfile, format="gtf", feature.type=NULL) 
setwd(path)
class(GTF)
attributes(GTF)

gene = "gene"
GTF1 = GTF[(elementMetadata(GTF)[,2] %in% gene)]
#set gene_id to names of GRRanges object
names(GTF1) <- elementMetadata(GTF1)[,7]
rm(GTF)

targets = read.table("results/consensus.clusters.txt",header = T,sep = "\t",na.strings =c("","NaN"),quote=NULL,stringsAsFactors=F,dec=".",fill=TRUE)
targets_rd = targets
#targets_rd$score <- 222
targets_rd <- targets_rd[,c(2:4,1,6,5)]
head(targets_rd)
targets_rd$consensus.cluster <- paste("cl_",as.character(targets_rd$consensus.cluster),sep="")
targets_rd <- BED2RangedData(targets_rd,header=F)
targets_rd = annotatePeakInBatch(targets_rd, AnnotationData = GTF1,select =  "first")
targets_rd <- as.data.frame(targets_rd)
row.names(targets_rd) <- sapply(strsplit(row.names(targets_rd),"\\."),"[[",1)
targets_rd$consensus.cluster = row.names(targets_rd)

targets$consensus.cluster = paste("cl_",as.character(targets$consensus.cluster),sep="")
annot <- merge(targets,targets_rd,by= "consensus.cluster")
write.table(annot,"results/consensus.clusters_annot.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,  col.names = NA, qmethod = c("escape", "double"))
write.table(targets,"results/consensus.clusters_new.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,  col.names = NA, qmethod = c("escape", "double"))

annot_sub = annot
kx = names(annot_sub) %in% c("consensus.cluster","feature")
annot_sub = annot_sub[kx]
names(annot_sub) = c("ID","gene_id")

r1 = read.table("results/consensus.clusters.normalized.tpm.txt",header = T,sep = "\t",na.strings =c("","NaN"),quote=NULL,stringsAsFactors=F,dec=".",fill=TRUE)
names(r1)
r1$ID = paste("cl_",r1$X,sep="")
r1$consensusclusterid = r1$X
r1$X <- NULL
r1 = r1[,c(ncol(r1)-1,ncol(r1),1:60)]
head(r1)
r2 = merge(r1,annot_sub,by="ID")
write.table(r2,"results/consensus.clusters.normalized.tpm.annot.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,  col.names = T, qmethod = c("escape", "double"))


#get plot about distribution across promoter region,introns,etc see cager vignette to plotannot()
#now plot distribution as in plotannot()
rm(list = ls(all = TRUE)); graphics.off()
path <- getwd(); setwd(path)
path

d1 = read.table("results/consensus.clusters_annot.txt",header = T,sep = "\t",na.strings =c("","NaN"),quote=NULL,stringsAsFactors=F,dec=".",fill=TRUE)
d1$ID = d1$consensus.cluster
#d1$consensus.cluster <- NULL
d1 = d1[,c("ID","insideFeature")]
d1$feature <- d1$insideFeature
d1$insideFeature <- NULL
d2 = read.table("results/consensus.clusters.normalized.tpm.txt",header = T,sep = "\t",na.strings =c("","NaN"),quote=NULL,stringsAsFactors=F,dec=".",fill=TRUE)
d2$ID = paste("cl_",d2$X,sep="")
d2$X <- NULL
me = merge(d1,d2,by="ID")
me = melt(me,id.vars = c("ID","feature"))
me$QBiC.Code =me$variable
m <- read.table("results/metadata.txt",header = T,sep = "\t",na.strings =c("","NaN"),quote=NULL,stringsAsFactors=F,dec=".",fill=TRUE)
me <- merge(me,m,by="QBiC.Code")
head(me)

pdf("feature_consensus_cluster_tpm_distribution.pdf")
ggplot(data=me, aes(x=QBiC.Code, y=value,fill=feature)) + geom_bar(position = "fill",stat = "identity") +
  ylab('relative contribution') +
  xlab('') +
  ggtitle('') +
  #facet_wrap(~chr) +
  coord_flip() +
  theme(text = element_text(size=6),
        axis.text.x = element_text(angle=45, vjust=1,hjust=1)) #vjust and hjust adjust the vertical and 
dev.off()


#hierarchical clustering
#clean up again
rm(list = ls(all = TRUE)); graphics.off()
path <- getwd(); setwd(path)
path

dat1 = read.table("results/consensus.clusters.normalized.tpm.txt",header = T,sep = "\t",na.strings =c("","NaN"),quote=NULL,stringsAsFactors=F,dec=".",fill=TRUE)
dat1$consensus.cluster = paste("cl_",dat1$X,sep="")
row.names(dat1) <- dat1$consensus.cluster
dat1$consensus.cluster <- NULL;dat1$X <- NULL

cl=(kmeans(dat1,4))
str(cl)
cl$size
cl$withinss

dat1$cluster=factor(cl$cluster)
centers=as.data.frame(cl$centers)
dat1$consensus.cluster = row.names(dat1)
dat1 <- melt(dat1,id.vars = c("cluster","consensus.cluster"))
dat1$QBiC.Code <- dat1$variable
m <- read.table("results/metadata.txt",header = T,sep = "\t",na.strings =c("","NaN"),quote=NULL,stringsAsFactors=F,dec=".",fill=TRUE)
dat1 <- merge(dat1,m,by="QBiC.Code")
head(dat1)

dat2 = data.frame(cluster = 1:4, size=cl$size)
dat1 = merge(dat1,dat2,by="cluster")
dat1$info = paste("cl",dat1$cluster,"size",dat1$size,sep=".")
head(dat1)
dat1$X <- paste(dat1$Condition..genotype,dat1$Condition..treatment,dat1$Condition..time,sep = "_")

pdf("consensus_clusters_kmeans_clustering_tpm_values1_log.pdf")
ggplot(dat1, aes(x=X, y=log2(value+0.1), fill=cluster)) +
  xlab('') +
  ylab('log(tpm)') +
  geom_violin(trim = F,draw_quantiles = c(0.25, 0.5, 0.75)) + 
  #coord_flip() +
  facet_wrap(~info,nrow = 1,labeller = label_parsed) +
  #scale_fill_manual(values=c("yellow", "darkgoldenrod","darkblue")) +
  theme(text = element_text(size=6),
        axis.text.x = element_text(angle=90, vjust=1,hjust=1))
dev.off()

system("mv *.pdf results/")
##system("mv *.bedGraph results/")
system("mv *.bed results/")
system("rm -rf results/bams")

###end part 1
#end of script
####-------------save Sessioninfo
fn <- paste("results/sessionInfo_",format(Sys.Date(), "%d_%m_%Y"),".txt",sep="")
sink(fn)
sessionInfo()
sink()
####---------END----save Sessioninfo

system("cp -r results_part6 results_all_backup")
system("mkdir tmp")
system("mv results_part* tmp/")
system("rm -rf tmp")
#end
