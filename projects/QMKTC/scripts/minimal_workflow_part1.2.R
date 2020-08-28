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

#check if genome is available
available.genomes()

##dir.create("results")
# 
# rm(ensembl)
# bo=0
# while(bo!=20){
#   ensembl = try(useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl"),silent=TRUE)
#   if (class(ensembl)=="try-error") {
#     cat("ERROR1: ", ensembl, "\n")
#     Sys.sleep(1)
#     print("reconnecting...")
#     bo <- bo+1
#     print(bo)
#   } else
#     break 
# }
# 
# save(list=c("ensembl"), file="ensembl.Rdata")
# rm(ensembl)


###################################################
### input files
###################################################
input <- list.files("results/bams",recursive = T, full.names = TRUE, pattern = ".out.bam")
input;basename(input)


#metadata
meta=0
meta = list.files(path,pattern = "sample_preparations")
meta = read.table(meta,header = T,sep = "\t",na.strings =c("","NaN"),quote=NULL,stringsAsFactors=F,dec=".",fill=TRUE)
#sort metadata by timepoint
meta$nn = input
meta$match <- NULL
f <- as.character(meta$QBiC.Code)
for (i in f) {
  print(i)
  meta[grepl(i, meta$nn), "match"] <- i
}
meta = subset(meta, ! match == 0 )
stopifnot(identical(meta$QBiC.Code,meta$match))

write.table(meta,"results/metadata.txt", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,  col.names = T, qmethod = c("escape", "double"))
#now ready to asssign a new vector with better filenames:
names_new <- as.character(meta$QBiC.Code)
#clean up meta
nk <- names(meta) %in% c("Extract.Code.s.","Extract.Name.s.","Source","Source.Name.s.","Source.Lab.ID.s.")
meta <- meta[!nk]


###################################################
### Creating a CAGEset object
###################################################
myCAGEset <- new("CAGEset", genomeName = "BSgenome.Hsapiens.UCSC.hg19",
                 inputFiles = input, inputFilesType = "bam",
                 sampleLabels = names_new)
sampleLabels(myCAGEset)
myCAGEset



###################################################
### now read data
###################################################
getCTSS(myCAGEset,removeFirstG = TRUE, correctSystematicG = TRUE)
#(Note that in case when a CAGEset object is created by coercion from a data.frame there is no need to call the above function to read in the data, as the data will be loaded into CAGEset during coercion)

stopifnot(identical(input,myCAGEset@inputFiles))
#TRUE
save.image("results/cage_data_1.RData")

system("cp -r results results_part2")
#end
