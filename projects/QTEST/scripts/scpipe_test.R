library(scPipe)
library(SingleCellExperiment)
data_dir = tempdir()

ERCCfa_fn = system.file("extdata", "ERCC92.fa", package = "scPipe")
ERCCanno_fn = system.file("extdata", "ERCC92_anno.gff3", package = "scPipe")
barcode_annotation_fn = system.file("extdata", "barcode_anno.csv", package = "scPipe")

fq_R1 = system.file("extdata", "simu_R1.fastq.gz", package = "scPipe")
fq_R2 = system.file("extdata", "simu_R2.fastq.gz", package = "scPipe")


sc_trim_barcode(file.path(data_dir, "combined.fastq.gz"),
                fq_R1,
                fq_R2,
                read_structure = list(bs1=-1, bl1=0, bs2=6, bl2=8, us=0, ul=6))
