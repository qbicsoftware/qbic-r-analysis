#clean up
rm(list = ls(all = TRUE)) # clear all variables
graphics.off()
#set path to location where script starts, for portability
path <- getwd()
setwd(path)
path


dir.create("results")
###################################################
### input files
###################################################
#filter input files
system("cp -r bams results/bams")
setwd("results/bams")
cmd="for i in *;do echo $i;samtools index $i;samtools view -o $i.out.bam $i chrX chrY `seq 1 22 | sed 's/^/chr/'`;done"
system(cmd)
setwd(path);path

system("cp -r results results_part1")
#end

