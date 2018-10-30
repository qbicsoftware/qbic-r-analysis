> sessionInfo()$otherPkgs
#!/usr/bin/env Rscript
#' ---
#' title: "Quality control analysis report of Microarray Data "
#' author: "author x"
#' ---
#SET THE WORKSPACE TO BE THE FOLDER OF THE PROJECT YOU ARE INTERESTED IN E.G. QMORAg
#source("https://bioconductor.org/biocLite.R")


if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)}
#'load packages needed for oligo run
library("siggenes")
library("RColorBrewer")
library("multtest")
library("limma")
library("oligo")
library("genefilter")
library("gplots")
library("ggplot2")
library("dendextend")
library("statmod")
library("annotate")

#+ clean up,echo=F
library(knitr)
knitr::opts_chunk$set(
  fig.width = 10,
  fig.height = 10,
  fig.path = 'qc_results/plots/',
  echo = T,
  warning = FALSE,
  message = FALSE
)

#+ clean_up,echo=F
rm(list = ls(all = TRUE)) #' clear all variables
graphics.off()
arguments = commandArgs(trailingOnly=TRUE)
path <- arguments[1]
<<<<<<< HEAD
#path = "/Users/Timo/Desktop/MITO_PD/QLSNK/qc_results"
=======

>>>>>>> qbicsoftware/master
setwd(path)

#' create directories needed
dir.create("qc_results")
dir.create("qc_results/raw_data")
dir.create("qc_results/plots")
dir.create("qc_results/tables")
dir.create("qc_results/metadata")
dir.create("qc_results/final")


if (file.exists("qc_results/metadata/metadata.txt")) {
  file.remove("qc_results/metadata/metadata.txt")
}


#' prepare raw and metadata files...
#+ prepare_files, echo=F



files = list.files(path, recursive = F, pattern = ".CEL")
for (i in files)
{
  cmd = paste("cp ", paste(i), " qc_results/raw_data/", sep = "")
  system(cmd)
}

setwd(path)

<<<<<<< HEAD
m <- list.files(path, recursive = F, pattern = "sample_preparation")
=======
m <- list.files(path, recursive = T, pattern = "sample_preparation")
>>>>>>> qbicsoftware/master
for (i in m)
{
  cmd = paste("cp ", paste(i), " qc_results/metadata/", sep = "")
  system(cmd)
}

m <- list.files("qc_results/metadata/")
<<<<<<< HEAD

=======
#RECALL: paste command adjusted for QCANG
>>>>>>> qbicsoftware/master
m <-
  read.table(
    paste("qc_results/metadata/", m, sep = ""),
    header = T,
    sep = "\t",
    na.strings = c("", "NaN"),
    quote = NULL,
    stringsAsFactors = F,
    dec = ".",
    fill = TRUE
  )
### create a new column `x` with all the columns collapsed together
cols <- names(m)
m$x <- apply(m[, cols], 1, paste, collapse = "_")
m$x <- gsub(" ", "", m$x)

#Rename cel files to fit the QBiC code obtained from the metadata for easy loading and copy them to cel_files folder
cel_names = m$QBiC.Code
celFiles <- list.celfiles("qc_results/raw_data", full.names = T)

if(!dir.exists("qc_results/raw_data/cel_files/")){
  dir.create("qc_results/raw_data/cel_files/")
}
for(i in 1:length(celFiles)){
  file_to_copy = grep(cel_names[i],celFiles,value=TRUE)
  copy_to = paste("qc_results/raw_data/cel_files/",cel_names[i],".CEL",sep="")
  file.copy(file_to_copy,copy_to)
}



###use filename to create a vector of regexes
files <-
  list.files("qc_results/raw_data/",
<<<<<<< HEAD
             recursive = F,
=======
             recursive = T,
>>>>>>> qbicsoftware/master
             pattern = ".CEL")
files <- gsub(".CEL", "", files)

#only take file identifier till 2nd underscore
files = paste(sapply(strsplit(files, "_"), "[[", 1), sapply(strsplit(files, "_"), "[[", 2), sep =
                "_")
m$xx <- 0
files

for (i in files) {
  ifelse(m[grepl(i, m$x), "xx"] <- i, 0)
}
m <- subset(m,!xx == 0)
m$xx <- paste(m$xx, ".CEL", sep = "")
m$filenames <- m$xx
m$xx <- NULL
m$x <- NULL

#' then write newly created metadata table

write.table(
  m,
  "qc_results/metadata/metadata.txt",
  append = FALSE,
  quote = FALSE,
  sep = "\t",
  eol = "\n",
  na = "NA",
  dec = ".",
  row.names = F,
  col.names = T,
  qmethod = c("escape", "double")
)

setwd(path)



#'Load metadata as AnnotatedDataFrame and load cel files
pd = read.AnnotatedDataFrame("qc_results/metadata/metadata.txt", header = TRUE)
row.names(pd) = pd[["filenames"]]
#'check filenames
row.names(pd)

#+ check,echo=F
setwd(path)

#' location of raw data cel files
celFiles <- list.celfiles("qc_results/raw_data/cel_files", full.names = T)
#stopifnot(identical(basename(celFiles), pd[["filenames"]]))

celFiles
num_cel <- length(celFiles)

#' then load them:
data <- read.celfiles(celFiles, phenoData = pd, verbose = T)

cel_counter <<- 1
cel_next <<- 15

#'a simple check:
identical(sampleNames(data), row.names(pd))


#' next fit Probe Level Models (PLMs) with probe- level and sample-level parameters.
#' The resulting object is an oligoPLM object, which stores parameter estimates, residuals and weights.

Pset <- fitProbeLevelModel(data)

#'Prepare color and groups for plots
#To find the column for labeling the script scans for a column that contains the tag control and uses this column as labels

conditions = grep("Condition", names(m), value = TRUE)

if (length(conditions) > 0) {
  d <- pData(data)[conditions]
  if (length(conditions) > 1) {
    grps <- paste(d[conditions][[1]], d[conditions][[2]], sep = "_")
  } else{
    grps <- d[conditions][[1]]
  }
  grps <- as.factor(grps)
  grps
  
  
}else{
  #WORKAROUND for QMTPD: Only do this for QMTPD project as conditions are not in separate column in this case
  expr = "PARKIN"
  cond = grep(expr,m$Source.Name.s.,value=TRUE)
  first = regexpr(expr,cond)
  keep1 = substr(cond,first,first+nchar(expr))
  
  expr = "LRRK2"
  cond = grep(expr,m$Source.Name.s.,value=TRUE)
  first = regexpr(expr,cond)
  keep2 = substr(cond,first,first+nchar(expr))
  
  expr = "CONTROL"
  cond = grep(expr,m$Source.Name.s.,value=TRUE)
  first = regexpr(expr,cond)
  keep3 = substr(cond,first,first+nchar(expr))
  grps = c(keep1,keep2,keep3)
  grps = as.factor(grps)
  #END WORKAROUND QMTPD
  if(length(grps) == 0){
  col = "white"
  grps = "sample"}
}

col = brewer.pal(length(levels(grps)), "Set1")
col = c(c(col)[grps])
#'now display groups and associated colours
grps
col


#' Basic quality analysis
#' 1) Correlation analysis between the samples and their repeated measurements
cor <- cor(exprs(data), use = "everything", method = c("pearson"))
write.table(
  cor,
  "qc_results/tables/pearson_correlation_all_data.tsv",
  append = FALSE,
  quote = FALSE,
  sep = "\t",
  eol = "\n",
  na = "NA",
  dec = ".",
  row.names = T,
  col.names = NA,
  qmethod = c("escape", "double")
)



#' 2) Boxplot of raw log intensities
ncol(exprs(data))
cel_counter <<- 1
cel_next <<- 15


svg("qc_results/plots/Boxplot_raw_intensities_%d.svg",onefile = FALSE)
for (i in 1:ceiling(length(celFiles) / 15)) {
  par(oma = c(3, 3, 3, 3), xpd = TRUE)
  #par(mfrow = c(3,3))
  
  if (col == "white") {
    boxplot(
      data[, cel_counter:cel_next],
      which = 'all',
      main = "Boxplot of the raw log intensities pre normalization",
      col = "white",
      xlab = "",
      main = "",
      ylab = "log2 signal intensity (PM+bg)",
      cex.axis = 0.5,
      las = 2
    )
  }
  else{
    boxplot(
      data[, cel_counter:cel_next],
      which = 'all',
      col = levels(factor(col)),
      xlab = "",
      main = "Boxplot of the raw log intensities pre normalization",
      ylab = "log2 signal intensity (PM+bg)",
      cex.axis = 0.5,
      las = 2
    )
    legend(
      "topright",
      col = levels(factor(col)),
      lwd = 1,
      cex = 0.5,
      legend = levels(grps),
      inset = c(0, 0)
    )
    
  }
  #plot.new()
  cel_counter <<- cel_counter + 15
  cel_next <<- cel_next + 15
  if(cel_next > length(celFiles)){
    cel_next <<- length(celFiles)
  }
}


dev.off()


#also make a boxplot using all samples at once for overview
svg("qc_results/plots/Boxplot_raw_intensities_all_samples.svg")
par(oma = c(3, 3, 3, 3), xpd = TRUE)
boxplot(
  data,
  which = 'all',
  col = levels(factor(col)),
  xlab = "",
  main = "Boxplot of the raw log intensities pre normalization",
  ylab = "log2 signal intensity (PM+bg)",
  cex.axis = 0.5,
  las = 2
)
legend(
  "topright",
  col = levels(factor(col)),
  lwd = 1,
  cex = 0.5,
  legend = levels(grps),
  inset = c(0, 0)
)

dev.off()

#' 3) Pseudo chip images
<<<<<<< HEAD
length = length(sampleNames(data))
for (chip in 1:length) {
  png(paste(
    "qc_results/plots/pseudo_image.",
    sampleNames(data)[chip],
    ".png",
    sep = ""
  ))
  image(data[, chip])
  dev.off()
}
=======
#length = length(sampleNames(data))

##
#for (chip in 1:length) {
#  png(paste(
#    "qc_results/plots/pseudo_image.",
#    sampleNames(data)[chip],
#    ".png",
#    sep = ""
#  ))
##  image(data[, chip])
#  dev.off()
#}
>>>>>>> qbicsoftware/master



##
# 4) RLE and NUSE plots

nuse.dat <- NUSE(Pset, type = "values")

cel_counter <- 1
cel_next <- 15

svg("qc_results/plots/NUSE_plot_%d.svg")
for (i in 1:ceiling(length(celFiles) / 15)) {
  par(oma = c(12, 3, 3, 3))
  boxplot(
    nuse.dat[,cel_counter:cel_next],
    main = "NUSE plot",
    ylim = c(0.5, 2),
    outline = FALSE,
    col = col,
    las = 2,
    cex.axis = 0.5,
    ylab = "Normalized Unscaled Error (NUSE) values",
    whisklty = "dashed",
    staplelty = 1,
    cex.axis = 0.75
  )
  abline(h=1.0,col="green")
  abline(h=1.1,col="red")
  
  cel_counter <<- cel_counter + 15
  cel_next <<- cel_next + 15
  if(cel_next > length(celFiles)){
    cel_next <<- length(celFiles)
  }
}
dev.off()

cel_counter <- 1
cel_next <- 15

rle.dat <- RLE(Pset, type = "values")

svg("qc_results/plots/RLE_plot_%d.svg")
for (i in 1:ceiling(length(celFiles) / 15)) {
  par(oma = c(12, 3, 3, 3))
  boxplot(
    rle.dat[,cel_counter:cel_next],
    main = "RLE plot",
    ylim = c(-8, 8),
    outline = FALSE,
    col = col,
    las = 2,
    cex.axis = 0.5,
    ylab = "Relative Log Expression (RLE) values",
    whisklty = "dashed",
    staplelty = 1,
    cex.axis = 0.75
  )
  cel_counter <<- cel_counter + 15
  cel_next <<- cel_next + 15
  if(cel_next > length(celFiles)){
    cel_next <<- length(celFiles)
  }
  }
dev.off()
##


#' 5) Histogram to compare log2 intensities vs density between arrays
#'

cel_counter <<- 1
cel_next <<- 15

#adjust color to be black for histograms
if(col == "white"){
  col = "black"
}

svg("qc_results/plots/Histogramm_log2_intensities_vs_density_%d.svg")
par(oma = c(12, 3, 3, 3))
par(mfrow = c(3, 3))



for (i in 1:ceiling(length(celFiles) / 15)) {
  hist(
    data[,cel_counter:cel_next],
    col = col,
    lty = 1,
    xlab = "log2 intensity",
    ylab = "density",
    main = "Histogramm of the log2 intensities vs density",
    xlim = c(2, 12),
    type = "l"
  )
  #legend("topright",col = col, lwd=1, legend=sampleNames(data),cex=0.5)
  cel_counter <<- cel_counter + 15
  cel_next <<- cel_next + 15
  if(cel_next > length(celFiles)){
    cel_next <<- length(celFiles)
  }
}
dev.off()

#make histogram for all samples
svg("qc_results/plots/Histogramm_log2_intensities_vs_density_all_samples.svg")
par(oma = c(12, 3, 3, 3))
par(mfrow = c(3, 3))
hist(
    data,
    col = col,
    lty = 1,
    xlab = "log2 intensity",
    ylab = "density",
    main = "Histogramm of the log2 intensities vs density",
    xlim = c(2, 12),
    type = "l"
  )
dev.off()



#If the color was set to be black before make it white again for the normalised boxplots
if(col == "black"){
  col = "white"
}

#' 6) MA plots raw data - only do MA plots for normalized data
#png("qc_results/plots/MA_plot_before_normalization_groups.png")
#par(oma=c(3,3,3,3))
#MAplot(data2, pairs=TRUE,na.rm=TRUE)
#dev.off()



#'7) PCA plot before normalization
#'
pca_before <-
  prcomp(t(exprs(data)),
         scores = TRUE,
         scale. = TRUE,
         cor = TRUE)
eigs <- pca_before$sdev^2
summary(pca_before)
proportions = eigs / sum(eigs)

#' sqrt of eigenvalues
pca_before$sdev
#'loadings
head(pca_before$rotation)
#'PCs (aka scores)
head(pca_before$x)

#' create data frame with scores
scores_before = as.data.frame(pca_before$x)
#' plot of observations

#'reorder grps just for pca plot
grps_pca <- grps

svg("qc_results/plots/PCA_before_normalization.svg")

if (col != "white") {
  ggplot(data = scores_before, aes(x = PC1, y = PC2, colour = grps_pca)) +
    geom_hline(yintercept = 0, colour = "gray65") +
    geom_vline(xintercept = 0, colour = "gray65") +
    #'geom_text(colour = "black",label=sampleNames(data), size = 2,angle=40) +
    #'scale_fill_manual(values=c("#'E41A1C", "#'377EB8", "#'4DAF4A"), breaks=c("Parabel", "Simbox", "Texus"), labels=c("Parabel", "Simbox", "Texus")) +
    geom_point(aes(shape = factor(grps)), size = 2) +
    scale_shape_manual(values = 1:nlevels(grps))  +
    #'scale_colour_manual(values = c("#'E41A1C","#'377EB8", "#'4DAF4A"))
    #'scale_shape_manual(values=1:nlevels(col)) +
    theme(legend.title = element_blank()) +  #'#' turn off legend title
    ggtitle("PCA plot before normalization")+
    xlab(paste("PC1 proportion of variance: ", round(proportions[1],digits=4)))+
    ylab(paste("PC2 proportion of variance: ", round(proportions[2],digits=4)))
}


if (col == "white") {
  ggplot(data = scores_before, aes(x = PC1, y = PC2, colour = "white")) +
    geom_hline(yintercept = 0, colour = "gray65") +
    geom_vline(xintercept = 0, colour = "gray65") +
    #'geom_text(colour = "black",label=sampleNames(data), size = 2,angle=40) +
    #'scale_fill_manual(values=c("#'E41A1C", "#'377EB8", "#'4DAF4A"), breaks=c("Parabel", "Simbox", "Texus"), labels=c("Parabel", "Simbox", "Texus")) +
    geom_point(aes(col = "white")) +
    #'scale_colour_manual(values = c("#'E41A1C","#'377EB8", "#'4DAF4A"))
    #'scale_shape_manual(values=1:nlevels(col)) +
    theme(legend.title = element_blank()) +  #'#' turn off legend title
    ggtitle("PCA plot before normalization")+
    xlab(paste("PC1 proportion of variance: ", round(proportions[1],digits=4)))+
    ylab(paste("PC2 proportion of variance: ", round(proportions[2],digits=4)))
  
}
dev.off()
 setwd(path)


#' Data Normalization of full data and subset for boxplots
#'
eset <- rma(data)

#'Quality plots after normalization
#'1) Boxplot after normalization
ncol(exprs(data))
cel_counter <<- 1
cel_next <<- 15

svg("qc_results/plots/Boxplot_after_normalization_%d.svg")
for (i in 1:ceiling(length(celFiles) / 15)) {
  par(oma = c(3, 3, 3, 3), xpd = TRUE)
  #par(mfrow = c(3,3))
  if (col == "white") {
    boxplot(
      eset[, cel_counter:cel_next],
      which = 'all',
      col = "white",
      xlab = "",
      main = "Boxplot of normalized log 2 intensities",
      ylab = "log2 signal intensity (PM+bg)",
      cex.axis = 0.5,
      las = 2
    )
  }
  else{
    boxplot(
      eset[,cel_counter:cel_next],
      which = 'all',
      col = levels(factor(col)),
      xlab = "",
      main = "Boxplot of normalized log 2 intensities",
      ylab = "log2 signal intensity (PM+bg)",
      cex.axis = 0.5,
      las = 2
    )
  legend(
    "topright",
    col = levels(factor(col)),
    lwd = 1,
    cex = 0.5,
    legend = levels(grps),
    inset = c(0, 0)
  )
  }
  cel_counter <<- cel_counter + 15
  cel_next <<- cel_next + 15
  if(cel_next > length(celFiles)){
    cel_next <<- length(celFiles)
  }
}

dev.off()

svg("qc_results/plots/Boxplot_after_normalization_all_samples.svg")
par(oma = c(3, 3, 3, 3), xpd = TRUE)

boxplot(
  eset,
  which = 'all',
  col = levels(factor(col)),
  xlab = "",
  main = "Boxplot of normalized log 2 intensities",
  ylab = "log2 signal intensity (PM+bg)",
  cex.axis = 0.5,
  las = 2
)
legend(
  "topright",
  col = levels(factor(col)),
  lwd = 1,
  cex = 0.5,
  legend = levels(grps),
  inset = c(0, 0)
)
dev.off()


#' 2) The MA plot also allows summarization, so groups can be compared more easily:
svg("qc_results/plots/MA_plot_after_normalization_groups.svg")
par(oma=c(3,3,3,3))
MAplot(exprs(eset),main="MA plot for each combination of groups" ,pairs = TRUE,groups=grps,cex=0.7)

dev.off()
 #' 3) Clustering
#'for arrays, problem: arrays are not row names but at column position, thus transpose is needed
d <- dist(t(exprs(eset))) #' find distance matrix
d
hc <- hclust(d)               #' apply hierarchical clustering

dend <- as.dendrogram(hc)
#'remember groups, assign a new color code
labels_colors(dend) <- col[order.dendrogram(dend)]
#'colorCodes = c("red", "blue", "green")
svg("qc_results/plots/Cluster_Dendogram.svg")
par(oma = c(10, 2, 2, 2))
dend %>% set("labels_cex", 0.5) %>% plot(main="Clustering of samples using normalized log 2 values")
legend(
  "topright",
  col = levels(factor(col)),
  lwd = 1,
  cex = 0.5,
  legend = levels(grps)
)
dev.off()


#' 4) PCA after normalization
pca <- prcomp(t(exprs(eset)), scores = TRUE, cor = TRUE)

summary(pca)
#' sqrt of eigenvalues
pca$sdev
#'loadings
head(pca$rotation)
#'PCs (aka scores)
head(pca$x)

eigs <- pca$sdev^2
summary(pca)
proportions = eigs / sum(eigs)


#' create data frame with scores
scores_after = as.data.frame(pca$x)
#' plot of observations
svg("qc_results/plots/PCA_after_normalization.svg")
if (col != "white") {
  ggplot(data = scores_after, aes(x = PC1, y = PC2, colour = grps)) +
    geom_hline(yintercept = 0, colour = "gray65") +
    geom_vline(xintercept = 0, colour = "gray65") +
    #'geom_text(colour = "black",label=sampleNames(data), size = 2,angle=40) +
    #'scale_fill_manual(values=c("#'E41A1C", "#'377EB8", "#'4DAF4A"), breaks=c("Parabel", "Simbox", "Texus"), labels=c("Parabel", "Simbox", "Texus")) +
    geom_point(aes(shape = factor(grps)), size = 2) +
    scale_shape_manual(values = 1:nlevels(grps))  +
    #'scale_colour_manual(values = c("#'E41A1C","#'377EB8", "#'4DAF4A"))
    #'scale_shape_manual(values=1:nlevels(col)) +
    theme(legend.title = element_blank()) +  #'#' turn off legend title
    ggtitle("PCA plot after normalization") +
    xlab(paste("PC1 proportion of variance: ", round(proportions[1],digits=4)))+
    ylab(paste("PC2 proportion of variance: ", round(proportions[2],digits=4)))
}
if (col == "white") {
  ggplot(data = scores_after, aes(x = PC1, y = PC2, colour = "white")) +
    geom_hline(yintercept = 0, colour = "gray65") +
    geom_vline(xintercept = 0, colour = "gray65") +
    #'geom_text(colour = "black",label=sampleNames(data), size = 2,angle=40) +
    #'scale_fill_manual(values=c("#'E41A1C", "#'377EB8", "#'4DAF4A"), breaks=c("Parabel", "Simbox", "Texus"), labels=c("Parabel", "Simbox", "Texus")) +
    geom_point(aes(col = "white")) +
    #'scale_colour_manual(values = c("#'E41A1C","#'377EB8", "#'4DAF4A"))
    #'scale_shape_manual(values=1:nlevels(col)) +
    theme(legend.title = element_blank()) +  #'#' turn off legend title
    ggtitle("PCA plot after normalization") +
    xlab(paste("PC1 proportion of variance: ", round(proportions[1],digits=4)))+
    ylab(paste("PC2 proportion of variance: ", round(proportions[2],digits=4)))
  
}
dev.off()


#biocLite("affycoretools")
#' 5) Scree plot to verify plotting of PC1 vs PC2
library("affycoretools")
svg("qc_results/plots/PCs.svg")
plotPCA(exprs(eset),
        main = "Principal component analysis (PCA)",
        screeplot = TRUE,
        outside = TRUE)
dev.off()

#'unload affy related packages again as analysis is focused on using oligo package function:
detach("package:affycoretools", unload = TRUE)
#'detach("package:affy", unload=TRUE)


#' 6) Non-specific filtering of data,could also be done with e.g. IQR
sds = rowSds(exprs(eset))
sh = shorth(sds)

svg("qc_results/plots/Histogram_of_sds.svg")
hist(sds, breaks=50, xlab="standard deviation")
abline(v=sh, col="blue", lwd=3, lty=2)
dev.off()
eset_filt_sds = eset[sds>=sh,]
# check dimensions of datasets before and after filtering
dim(exprs(eset))
dim(exprs(eset_filt_sds))

#'write to file for customer to plot e.g. in excel:
write.exprs(eset, "qc_results/tables/RMAnorm_nonfiltered.txt")
#write.exprs(eset_filt_sds, "qc_results/tables/RMAnorm_sds.filtered.txt")


#'save esets
save(list = c("pd", "eset"), file = "qc_results/final/eset.Rdata")

#'save image
setwd(path)



#' end of script
#' save Sessioninfo
fn <-
  paste("qc_results/tables/sessionInfo_",
        format(Sys.Date(), "%d_%m_%Y"),
        ".txt",
        sep = "")
sink(fn)
sessionInfo()
sink()