#!/usr/bin/env Rscript

library(alakazam)
library(ggplot2)
library(data.table)
library(dplyr)
library(tigger)
library(shazam)
library(igraph)
library(gplots)
library(stringr)

theme_set(theme_bw(base_family = "ArialMT") + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(family="ArialMT")))



#datadir <- "Data"
datadir <- "../germlines"
outdir <- "repertoire_comparison"

# setwd to results folder (containing alakazam, shazam, etc. folders)

### Read all the tables as produced by the pipeline in the current folder and joins them together in the df_all dataframe

all_files <- system(paste0("find '",datadir,"' -name '*.tab'"), intern=T)

dir.create(outdir)
diversity_dir <- paste(outdir, "Diversity", sep="/")
abundance_dir <- paste(outdir, "Abundance", sep="/")
isotype_dir <- paste(outdir, "Isotype", sep="/")
vfamily_dir <- paste(outdir, "V_family", sep="/")
mutation_dir <- paste(outdir, "Mutational_load", sep="/")
dir.create(diversity_dir)
dir.create(abundance_dir)
dir.create(isotype_dir)
dir.create(vfamily_dir)
dir.create(mutation_dir)

# Set number of bootrstraps
nboot = 1000

# Generate one big dataframe from all patient dataframes
df_all = data.frame()
for (file in all_files){
   fname = file
   print(fname)

   df_pat <- read.csv(fname, sep="\t")
   
   df_all <- rbind(df_all, df_pat)

}
write.table(df_all, "all_data.tab", sep = "\t", quote=F, row.names = F, col.names = T)

#df_all <- read.csv("all_data.tab", sep = "\t")

# Remove underscores in these columns
df_all$TREATMENT <- sapply(df_all$TREATMENT, function(x) str_replace(as.character(x), "_", ""))
df_all$SOURCE <- sapply(df_all$SOURCE, function(x) str_replace(as.character(x), "_", ""))
df_all$EXTRACT_TIME <- sapply(df_all$EXTRACT_TIME, function(x) str_replace(as.character(x), "_", ""))
df_all$POPULATION <- sapply(df_all$POPULATION, function(x) str_replace(as.character(x), "_", ""))

# Annotate sample and samplepop (sample + population) by add ing all the conditions
df_all$SAMPLE <- as.factor(paste(df_all$TREATMENT, df_all$EXTRACT_TIME, df_all$SOURCE, sep="_"))
df_all$SAMPLE_POP <- as.factor(paste(df_all$TREATMENT, df_all$EXTRACT_TIME, df_all$SOURCE, df_all$POPULATION, sep="_"))

##########################
## ABUNDANCE PER PATIENT
#########################

# Per patient
abund <- estimateAbundance(df_all, clone="CLONE", group = "SAMPLE", ci=0.95, nboot=nboot)
abund@abundance$TREATMENT <- sapply(abund@abundance$SAMPLE, function(x) unlist(strsplit(as.character(x), "_"))[1])
abund@abundance$TIME_POINT <- sapply(abund@abundance$SAMPLE, function(x) unlist(strsplit(as.character(x), "_"))[2])
abund@abundance$PATIENT <- sapply(abund@abundance$SAMPLE, function(x) unlist(strsplit(as.character(x), "_"))[3])

abund_main <- paste0("Clonal abundance (n=", abund@n[1], ")")

p1 <- ggplot(abund@abundance, aes(x = rank, y = p, 
                             group = SAMPLE)) + 
  geom_ribbon(aes(ymin = lower, 
                  ymax = upper, fill = PATIENT), alpha = 0.4) + 
  geom_line(aes(color = PATIENT)) +
  ggtitle(abund_main) + 
  xlab("Rank") + ylab("Abundance") + 
  scale_x_log10(limits = NULL, 
                breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_continuous(labels = scales::percent) +
  facet_grid(cols = vars(TREATMENT), scales="free", drop = T)
ggsave(plot=p1, filename = paste0(abundance_dir,"/Rank_abundance_patient.svg"), device="svg", width = 25, height = 10, units="cm")

write.table(abund@abundance, file = paste0(abundance_dir, "/Abundance_data_patient.tsv"), sep="\t", quote = F, row.names = F)

########################
## DIVERSITY PER PATIENT
#######################

print("Diversity calculation")
# Plotting sample diversity per patient
sample_div <- alphaDiversity(abund, group="SAMPLE", min_q=0, max_q=4, step_q=0.05,
                             ci=0.95, nboot=nboot)
sample_main <- paste0("Sample diversity (n=", sample_div@n[1], ")")

sample_div@diversity$TREATMENT <- sapply(sample_div@diversity$SAMPLE, function(x) unlist(strsplit(as.character(x), "_"))[1])
sample_div@diversity$TIME_POINT <- sapply(sample_div@diversity$SAMPLE, function(x) unlist(strsplit(as.character(x), "_"))[2])
sample_div@diversity$PATIENT <- sapply(sample_div@diversity$SAMPLE, function(x) unlist(strsplit(as.character(x), "_"))[3])

p2 <- ggplot(sample_div@diversity, aes(x = q, y = d, 
                                       group = SAMPLE)) + 
  geom_ribbon(aes(ymin = d_lower, 
                  ymax = d_upper, 
                  fill = PATIENT), alpha = 0.4) +
  geom_line(aes(color = PATIENT)) +
  xlab("q") + ylab(expression(""^q * D)) +
  ggtitle(sample_main) + 
  facet_grid(cols=vars(TREATMENT))
ggsave(paste0(diversity_dir,"/Diversity_patient_grid.svg"), device="svg", width = 25, height = 10, units="cm")
ggsave(paste0(diversity_dir,"/Diversity_patient_grid.pdf"), device="pdf", width = 25, height = 10, units="cm")

# DIVERSITY AT Q=1

sample_div_q1 <- sample_div@diversity[which(sample_div@diversity$q == 1),]

sample_main <- paste0("Sample diversity at Q=1 (n=", sample_div@n[1], ")")

dodge <- position_dodge(width = 0.9)
g1 <- ggplot(sample_div_q1, aes(y=d, x=PATIENT)) + 
  geom_point(position=dodge, stat="identity") +
  geom_errorbar(aes(ymin=d-d_sd, ymax=d+d_sd), width = .2, position=dodge) +
  xlab("") + ylab("Diversity (q=1)") +
  ggtitle(sample_main) +
  facet_grid(cols=vars(TREATMENT), drop=T, space="free", scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(plot = g1, filename = paste0(diversity_dir,"/Diversity_Q1_test_patient.svg"), device="svg", 
       width = 20, height = 10, units="cm")
ggsave(plot = g1, filename = paste0(diversity_dir,"/Diversity_Q1_test_patient.pdf"), device="pdf", 
       width = 20, height = 10, units="cm")

# DIVERSITY AT Q=0

sample_div_q0 <- sample_div@diversity[which(sample_div@diversity$q == 0),]

sample_main <- paste0("Sample diversity at Q=0 (n=", sample_div@n[1], ")")

dodge <- position_dodge(width = 0.9)
g1 <- ggplot(sample_div_q1, aes(y=d, x=PATIENT)) + 
  geom_point(position=dodge, stat="identity") +
  geom_errorbar(aes(ymin=d-d_sd, ymax=d+d_sd), width = .2, position=dodge) +
  xlab("") + ylab("Diversity (q=0)") +
  ggtitle(sample_main) +
  facet_grid(cols=vars(TREATMENT), drop=T, space="free", scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(plot = g1, filename = paste0(diversity_dir,"/Diversity_Q0_test_patient.svg"), device="svg", 
       width = 20, height = 10, units="cm")
ggsave(plot = g1, filename = paste0(diversity_dir,"/Diversity_Q0_test_patient.pdf"), device="pdf", 
       width = 20, height = 10, units="cm")

sample_test <- sample_div@tests
write.table(sample_test, file = paste0(diversity_dir, "/Diversity_tests_data_patient.tsv"), 
            sep="\t", quote = F, row.names = F)

##########################
# ABUNDANCE PER POPULATION
##########################


abund <- estimateAbundance(df_all, clone="CLONE", group = "SAMPLE_POP", ci=0.95, nboot=nboot)
abund@abundance$TREATMENT <- sapply(abund@abundance$SAMPLE_POP, function(x) unlist(strsplit(as.character(x), "_"))[1])
abund@abundance$TIME_POINT <- sapply(abund@abundance$SAMPLE_POP, function(x) unlist(strsplit(as.character(x), "_"))[2])
abund@abundance$PATIENT <- sapply(abund@abundance$SAMPLE_POP, function(x) unlist(strsplit(as.character(x), "_"))[3])
abund@abundance$POPULATION <- sapply(abund@abundance$SAMPLE_POP, function(x) unlist(strsplit(as.character(x), "_"))[4])

abund_main <- paste0("Clonal abundance (n=", abund@n[1], ")")
p1 <- ggplot(abund@abundance, aes(x = rank, y = p, 
                             group = SAMPLE_POP)) + 
  geom_ribbon(aes(ymin = lower, 
                  ymax = upper, fill = POPULATION), alpha = 0.4) + 
  geom_line(aes(color = POPULATION)) +
  ggtitle(abund_main) + 
  xlab("Rank") + ylab("Abundance") + 
  scale_x_log10(limits = NULL, 
                breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(vars(PATIENT), scales="free", drop = T)
ggsave(plot=p1, filename = paste0(abundance_dir,"/Rank_abundance_patient_population.svg"), device="svg", 
       width = 30, height = 20, units="cm")

write.table(abund@abundance, file = paste0(abundance_dir, "/Abundance_data_patient_population.tsv"), sep="\t", 
            quote = F, row.names = F)

##########################
# ABUNDANCE PER STATE
##########################

abund_main <- paste0("Clonal abundance (n=", abund@n[1], ")")
p1 <- ggplot(abund@abundance, aes(x = rank, y = p, 
                                  group = SAMPLE_POP)) + 
  geom_ribbon(aes(ymin = lower, 
                  ymax = upper, fill = PATIENT), alpha = 0.4) + 
  geom_line(aes(color = PATIENT)) +
  ggtitle(abund_main) + 
  xlab("Rank") + ylab("Abundance") + 
  scale_x_log10(limits = NULL, 
                breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_continuous(labels = scales::percent) +
  facet_grid(rows=vars(POPULATION), cols=vars(TREATMENT), drop=T, scales="free") +
  #facet_wrap(vars(POPULATION,TREATMENT), drop = T, nrow=4, ncol=4) + 
  theme(strip.text.y = element_text(size = 5))
ggsave(plot=p1, filename = paste0(abundance_dir,"/Rank_abundance_population_state.svg"), device="svg", 
       width = 15, height = 20, units="cm")
ggsave(plot=p1, filename = paste0(abundance_dir,"/Rank_abundance_population_state.png"), device="png",
       width = 15, height = 20, units="cm")

###########################
# DIVERSITY PER POPULATION
###########################

# Plotting sample diversity for all populations
print("Diversity calculation population")

sample_div <- alphaDiversity(abund, group="SAMPLE_POP", min_q=0, max_q=4, step_q=0.05,
                             ci=0.95, nboot=nboot)
sample_main <- paste0("Sample diversity (n=", sample_div@n[1], ")")

sample_div@diversity$TREATMENT <- sapply(sample_div@diversity$SAMPLE_POP, function(x) unlist(strsplit(as.character(x), "_"))[1])
sample_div@diversity$TIME_POINT <- sapply(sample_div@diversity$SAMPLE_POP, function(x) unlist(strsplit(as.character(x), "_"))[2])
sample_div@diversity$PATIENT <- sapply(sample_div@diversity$SAMPLE_POP, function(x) unlist(strsplit(as.character(x), "_"))[3])
sample_div@diversity$POPULATION <- sapply(sample_div@diversity$SAMPLE_POP, function(x) unlist(strsplit(as.character(x), "_"))[4])

p2 <- ggplot(sample_div@diversity, aes(x = q, y = d, 
                                   group = SAMPLE_POP)) + 
  geom_ribbon(aes(ymin = d_lower, 
            ymax = d_upper, fill = POPULATION), alpha = 0.4) +
  geom_line(aes(color = POPULATION)) +
  xlab("q") + ylab(expression(""^q * D)) +
  ggtitle(sample_main) + 
  facet_wrap(vars(PATIENT), scales = "free_x")
ggsave(plot = p2, filename = paste0(diversity_dir,"/Diversity_patient_population.svg"), device="svg", 
width = 25, height = 20, units="cm")
ggsave(plot = p2, filename = paste0(diversity_dir,"/Diversity_patient_population.pdf"), device="pdf", 
width = 25, height = 20, units="cm")

###########################
# DIVERSITY PER STATE
###########################
p2 <- ggplot(sample_div@diversity, aes(x = q, y = d, 
                                       group = SAMPLE_POP)) + 
  geom_ribbon(aes(ymin = d_lower, 
                  ymax = d_upper, fill = PATIENT), alpha = 0.4) +
  geom_line(aes(color = PATIENT)) +
  xlab("q") + ylab(expression(""^q * D)) +
  ggtitle(sample_main) +
  theme(strip.text.y = element_text(size = 5)) +   
  facet_grid(rows=vars(POPULATION), cols=vars(TREATMENT), drop=T, scales="free_x")
#facet_wrap(vars(POPULATION,TREATMENT), scales = "free_x")
ggsave(plot = p2, filename = paste0(diversity_dir,"/Diversity_population_state.svg"), device="svg", 
       width = 30, height = 20, units="cm")
ggsave(plot = p2, filename = paste0(diversity_dir,"/Diversity_population_state.pdf"), device="pdf", 
       width = 30, height = 20, units="cm")

# Tests sample diversity for significance
print("Diversity calculation tests population")

# DIVERSITY AT Q=1

sample_div_q1 <- sample_div@diversity[which(sample_div@diversity$q == 1),]

paste0("Sample diversity at Q=1 (n=", sample_div@n[1], ")")

dodge <- position_dodge(width = 0.9)
g1 <- ggplot(sample_div_q1, aes(y=d, x=POPULATION, fill=POPULATION)) + 
  geom_point(position=dodge, stat="identity", aes(fill=POPULATION)) +
  geom_errorbar(aes(ymin=d-d_sd, ymax=d+d_sd), width = .2, position=dodge) +
  xlab("") + ylab("Diversity (q=1)") +
  ggtitle(sample_main) +
  facet_grid(cols=vars(PATIENT)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(plot = g1, filename = paste0(diversity_dir,"/Diversity_Q1_test_patient_population.svg"), device="svg", 
       width = 30, height = 10, units="cm")
ggsave(plot = g1, filename = paste0(diversity_dir,"/Diversity_Q1_test_patient_population.pdf"), device="pdf", 
       width = 30, height = 10, units="cm")

#Diversity Q1 by state
dodge <- position_dodge(width = 0.9)
g2 <- ggplot(sample_div_q1, aes(y=d, x=TREATMENT)) + 
  geom_jitter(width=0.05,stat="identity", aes(color=PATIENT)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.3, size = 0.3) +
  xlab("") + ylab("Diversity (q=1)") +
  ggtitle(sample_main) +
  theme(axis.text.x = element_text(angle = 49, hjust = 1), legend.position = "right")  +
  facet_wrap(vars(POPULATION), scales="free", drop = T)

ggsave(plot = g2, filename = paste0(diversity_dir,"/Diversity_Q1_population_state.svg"), device="svg", 
       width = 30, height = 20, units="cm")
ggsave(plot = g2, filename = paste0(diversity_dir,"/Diversity_Q1_population_state.pdf"), device="pdf", 
       width = 30, height = 20, units="cm")

# DIVERSITY AT Q=0

sample_div_q0 <- sample_div@diversity[which(sample_div@diversity$q == 0),]

paste0("Sample diversity at Q=0 (n=", sample_div@n[1], ")")

dodge <- position_dodge(width = 0.9)
g1 <- ggplot(sample_div_q0, aes(y=d, x=POPULATION, fill=POPULATION)) + 
  geom_point(position=dodge, stat="identity", aes(fill=POPULATION)) +
  geom_errorbar(aes(ymin=d-d_sd, ymax=d+d_sd), width = .2, position=dodge) +
  xlab("") + ylab("Diversity (q=0)") +
  ggtitle(sample_main) +
  facet_grid(cols=vars(PATIENT)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(plot = g1, filename = paste0(diversity_dir,"/Diversity_Q0_test_patient_population.svg"), device="svg", 
       width = 30, height = 10, units="cm")
ggsave(plot = g1, filename = paste0(diversity_dir,"/Diversity_Q0_test_patient_population.pdf"), device="pdf", 
       width = 30, height = 10, units="cm")

dodge <- position_dodge(width = 0.9)
g2 <- ggplot(sample_div_q0, aes(y=d, x=TREATMENT)) + 
  geom_jitter(width=0.05, stat="identity", aes(color=PATIENT)) +
  # geom_errorbar(aes(ymin=d-d_sd, ymax=d+d_sd), width = .2, position=dodge) +
  xlab("") + ylab("Diversity (q=0)") +
  ggtitle(sample_main) +
  facet_wrap(vars(POPULATION), scales="free_x", drop = T) +
  # facet_grid(cols=vars(TREATMENT), rows=vars(POPULATION)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.3, size = 0.3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right") 

ggsave(plot = g2, filename = paste0(diversity_dir,"/Diversity_Q0_population_state.svg"), device="svg", 
       width = 30, height = 15, units="cm")
ggsave(plot = g2, filename = paste0(diversity_dir,"/Diversity_Q0_population_state.pdf"), device="pdf", 
       width = 30, height = 15, units="cm")


sample_test <- sample_div@tests
write.table(sample_test, file = paste0(diversity_dir, "/Diversity_tests_data_patient_population.tsv"), 
            sep="\t", quote = F, row.names = F)


############
## ISOTYPES
############

print("Isotypes calculation")
# Plotting Isotype percentages per patient
df_all$ISOTYPE <- df_all$C_PRIMER

res <- df_all %>% group_by(ISOTYPE,SAMPLE,SOURCE,TREATMENT,EXTRACT_TIME) %>% dplyr::summarise(Seqs_isotype=n())
res <- with(res, res[order(SOURCE),])
res_sample <- df_all %>% group_by(SAMPLE) %>% dplyr::summarise(Seqs_total=n())

freqs <- merge(x=res, y=res_sample, all.x = T, by.x = "SAMPLE", by.y = "SAMPLE")
freqs$Freq <- (freqs$Seqs_isotype/freqs$Seqs_total)

g4 <- ggplot(freqs, aes(fill=ISOTYPE, y=Freq, x=ISOTYPE)) +
  geom_bar(position = "dodge", stat="identity") +
  xlab("") + 
  ylab("Frequency") +
  ggtitle("Isotype frequency") +
  facet_wrap(vars(SOURCE), scales="free_x", drop = T) +
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))
ggsave(plot=g4, filename = paste0(isotype_dir,"/Isotype_frequencies_patient.svg"), device = "svg", 
  width = 18, height = 15, units = "cm")
ggsave(plot=g4, filename = paste0(isotype_dir,"/Isotype_frequencies_patient.pdf"), device = "pdf", 
  width = 18, height = 15, units = "cm")

write.table(freqs, file = paste0(isotype_dir,"/Isotype_frequencies_data.tsv"), sep="\t", quote=F, row.names = F)

# Plotting Isotype percentages per population
print("Isotypes calculation per population")
res <- df_all %>% group_by(ISOTYPE,SAMPLE_POP,SOURCE,TREATMENT,EXTRACT_TIME, POPULATION) %>% dplyr::summarise(Seqs_isotype=n())
res <- with(res, res[order(SOURCE),])
res_sample <- df_all %>% group_by(SAMPLE_POP) %>% dplyr::summarise(Seqs_total=n())

freqs <- merge(x=res, y=res_sample, all.x = T, by.x = "SAMPLE_POP", by.y = "SAMPLE_POP")
freqs$Freq <- (freqs$Seqs_isotype/freqs$Seqs_total)

g4 <- ggplot(freqs, aes(fill=ISOTYPE, y=Freq, x=ISOTYPE)) +
 geom_bar(position = "dodge", stat="identity") +
 xlab("") + 
 ylab("Frequency") +
 ggtitle("Isotype frequency") +
 facet_grid(cols=vars(SOURCE), rows=vars(POPULATION)) +
 theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))
ggsave(g4, filename = paste0(isotype_dir,"/Isotype_percentages_population.svg"), device = "svg", 
  width = 25, height = 20, units = "cm")
ggsave(g4, filename = paste0(isotype_dir,"/Isotype_percentages_population.pdf"), device = "pdf", 
  width = 25, height = 20, units = "cm")

write.table(freqs, file = paste0(isotype_dir, "/Isotype_frequencies_population_data.tsv"), sep="\t", quote = F, row.names = F)

#Plotting Isotype percentages per state
g4 <- ggplot(freqs, aes(y=Freq, x=TREATMENT)) +
  geom_point(position = "dodge", stat="identity", aes(color=SOURCE)) +
  xlab("") + 
  ylab("Frequency") +
  ggtitle("Isotype frequency") +
  stat_summary(fun = median, geom = "crossbar", width = 0.6, size = 0.3) +
  facet_grid(cols =vars(ISOTYPE), rows=vars(POPULATION),scales="free", drop = T)+
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))
ggsave(g4, filename = paste0(isotype_dir,"/Isotype_percentages_population_state.svg"), device = "svg", 
       width = 25, height = 12, units = "cm")
ggsave(g4, filename = paste0(isotype_dir,"/Isotype_percentages_population_state.pdf"), device = "pdf", 
       width = 25, height = 12, units = "cm")

##################
# V FAMILY USAGE              
##################

# Quantify V family usage by patient
print("V family usage calculation")
family <- countGenes(df_all, gene="V_CALL", groups="SAMPLE", 
                     mode="family")
family$TREATMENT <- sapply(family$SAMPLE, function(x) unlist(strsplit(as.character(x), "_"))[1])
family$TIME_POINT <- sapply(family$SAMPLE, function(x) unlist(strsplit(as.character(x), "_"))[2])
family$PATIENT <- sapply(family$SAMPLE, function(x) unlist(strsplit(as.character(x), "_"))[3])

g2 <- ggplot(family, aes(x=gene, y=seq_freq)) +
  scale_color_brewer(palette="Set1") +
  geom_point(aes(color=gene), size=4, alpha=0.8) +
  ggtitle("V Gene Family Usage") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  ylab("Frequency") +
  xlab("") +
  facet_wrap(vars(PATIENT), scales="free_x")
ggsave(filename = paste0(vfamily_dir, "/V_Family_distribution_patient.svg"), plot = g2, width = 18, height = 15, units = "cm")
ggsave(filename = paste0(vfamily_dir, "/V_Family_distribution_patient.pdf"), plot = g2, width = 18, height = 15, units = "cm")

write.table(family, file = paste0(vfamily_dir, "/V_family_distribution_data.tsv"), sep = "\t", quote = F, row.names = F)

# Quantify V family usage by population
print("V family usage calculation population")
family <- countGenes(df_all, gene="V_CALL", groups="SAMPLE_POP", 
                     mode="family")
family$TREATMENT <- sapply(family$SAMPLE_POP, function(x) unlist(strsplit(as.character(x), "_"))[1])
family$TIME_POINT <- sapply(family$SAMPLE_POP, function(x) unlist(strsplit(as.character(x), "_"))[2])
family$PATIENT <- sapply(family$SAMPLE_POP, function(x) unlist(strsplit(as.character(x), "_"))[3])
family$POPULATION <- sapply(family$SAMPLE_POP, function(x) unlist(strsplit(as.character(x), "_"))[4])

g2 <- ggplot(family, aes(x=gene, y=seq_freq)) +
  scale_color_brewer(palette="Set1") +
  geom_point(aes(color=gene), size=3, alpha=0.8) +
  ggtitle("V gene Family Usage") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  ylab("Frequency") +
  xlab("") +
  facet_grid(cols=vars(PATIENT), rows=vars(POPULATION))
ggsave(filename = paste0(vfamily_dir,"/V_Family_distribution_patient_population.svg"), plot = g2, 
  width = 30, height = 20, units = "cm")
ggsave(filename = paste0(vfamily_dir,"/V_Family_distribution_patient_population.pdf"), plot = g2, 
  width = 30, height = 20, units = "cm")

write.table(family, file = paste0(vfamily_dir, "/V_family_distribution_data_population.tsv"), sep = "\t", 
  quote = F, row.names = F)

# Quantify V family usage by state
g2 <- ggplot(family, aes(x=TREATMENT, y=seq_freq)) +
  scale_color_brewer(palette="Set1") +
  geom_point(aes(color=PATIENT), size=3, alpha=0.8) +
  ggtitle("V gene Family Usage") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  ylab("Frequency") +
  xlab("") +
  facet_grid(cols=vars(POPULATION), rows=vars(gene),scales="free", drop = T) 
#facet_grid(cols=vars(PATIENT,TREATMENT), rows=vars(POPULATION))
ggsave(filename = paste0(vfamily_dir,"/V_Family_distribution_patient_population_state.svg"), plot = g2, 
       width = 30, height = 20, units = "cm")
ggsave(filename = paste0(vfamily_dir,"/V_Family_distribution_patient_population_state.pdf"), plot = g2, 
       width = 30, height = 20, units = "cm")


#######################
# MUTATIONAL FREQUENCY
#######################

# Quantify mutational load by patient
print("Mutational load calculation")
# Calculate mutation counts
df_all_mut_counts <- observedMutations(df_all, sequenceColumn="SEQUENCE_IMGT",
                                     germlineColumn="GERMLINE_IMGT_D_MASK",
                                     regionDefinition=NULL,
                                     frequency=FALSE,
                                     combine=TRUE,
                                     nproc=4) 
# Calculate mutation frequencies
df_all_mut_freq <- observedMutations(df_all_mut_counts, sequenceColumn = "SEQUENCE_IMGT",
                                     germlineColumn = "GERMLINE_IMGT_D_MASK",
                                     regionDefinition = NULL,
                                     frequency = TRUE,
                                     combine=TRUE,
                                     nproc = 4)

mut_counts_freqs <- select(df_all_mut_freq, SEQUENCE_ID, SOURCE, EXTRACT_TIME, TREATMENT, POPULATION, SAMPLE, SAMPLE_POP, MU_COUNT, MU_FREQ)

res_mut <- mut_counts_freqs %>% group_by(SAMPLE,SOURCE,TREATMENT,EXTRACT_TIME) %>% 
  dplyr::summarise(MUTATION_MEAN_COUNT=mean(MU_COUNT), 
  MUTATION_MEDIAN_COUNT=median(MU_COUNT), MUTATION_SD_COUNT=sd(MU_COUNT), MUTATION_MEAN_FREQ=mean(MU_FREQ), MUTATION_MEDIAN_FREQ=median(MU_FREQ), MUTATION_SD_FREQ=mean(MU_FREQ), N_SEQS = n())
write.table(res_mut, file = paste0(mutation_dir,"/Mutation_stats_patient.tsv"), sep="\t", col.names = T, row.names = F, quote = F)

res_mut_pop <- mut_counts_freqs %>% group_by(SAMPLE_POP,SOURCE,TREATMENT,EXTRACT_TIME, POPULATION) %>% 
  dplyr::summarise(MUTATION_MEAN_COUNT=mean(MU_COUNT), MUTATION_MEDIAN_COUNT=median(MU_COUNT), MUTATION_SD_COUNT=sd(MU_COUNT), 
  MUTATION_MEAN_FREQ=mean(MU_FREQ), MUTATION_MEDIAN_FREQ=median(MU_FREQ), MUTATION_SD_FREQ=mean(MU_FREQ), N_SEQS = n())
write.table(res_mut_pop, file = paste0(mutation_dir,"/Mutation_stats_patient_population.tsv"), sep="\t", col.names = T, row.names = F, quote = F)


plot_mut_num <- ggplot(mut_counts_freqs, aes(fill=TREATMENT, y=MU_COUNT, x=SOURCE)) +
  geom_boxplot() +
  xlab("") + 
  ylab("Mutation Counts") +
  ggtitle("Mutation Counts") +
  facet_grid(cols=vars(TREATMENT), scales = "free", drop = T) +
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))
ggsave(plot=plot_mut_num, filename = paste0(mutation_dir,"/Mutation_count_patient.svg"), device = "svg", width = 25, height = 10, units = "cm")
ggsave(plot=plot_mut_num, filename = paste0(mutation_dir,"/Mutation_count_patient.png"), device = "png", width = 25, height = 10, units = "cm")


plot_mut_freq <- ggplot(mut_counts_freqs, aes(fill=TREATMENT, y=MU_FREQ, x=SOURCE)) +
  geom_boxplot() +
  xlab("") + 
  ylab("Mutation Frequency") +
  ggtitle("Mutation Frequency") +
  facet_grid(cols=vars(TREATMENT), scales = "free", drop = T) +
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))
ggsave(plot=plot_mut_freq, filename = paste0(mutation_dir,"/Mutation_frequency_patient.svg"), device = "svg", width = 25, height = 10, units = "cm")
ggsave(plot=plot_mut_freq, filename = paste0(mutation_dir,"/Mutation_frequency_patient.png"), device = "png", width = 25, height = 10, units = "cm")


plot_mut_num <- ggplot(mut_counts_freqs, aes(fill=POPULATION, y=MU_COUNT, x=POPULATION)) +
  geom_boxplot() +
  xlab("") + 
  ylab("Mutation Counts") +
  ggtitle("Mutation Counts per Population") +
  facet_wrap(vars(SOURCE), scales = "free_x", drop = T) +
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))
ggsave(plot=plot_mut_num, filename = paste0(mutation_dir,"/Mutation_count_patient_population.svg"), device = "svg", width = 25, height = 20, units = "cm")
ggsave(plot=plot_mut_num, filename = paste0(mutation_dir,"/Mutation_count_patient_population.png"), device = "png", width = 25, height = 20, units = "cm")


plot_mut_freq <- ggplot(mut_counts_freqs, aes(fill=POPULATION, y=MU_FREQ, x=POPULATION)) +
  geom_boxplot() +
  xlab("") + 
  ylab("Mutation Frequency") +
  ggtitle("Mutation Frequency per Population") +
  facet_wrap(vars(SOURCE), scales = "free_x", drop = T) +
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))
ggsave(plot=plot_mut_freq, filename = paste0(mutation_dir,"/Mutation_frequency_patient_population.svg"), device = "svg", width = 25, height = 20, units = "cm")
ggsave(plot=plot_mut_freq, filename = paste0(mutation_dir,"/Mutation_frequency_patient_population.png"), device = "png", width = 25, height = 20, units = "cm")

n_fun <- function(x){
  return(data.frame(y = 1.7*70,
                    label = length(x)))
}

# Calculate mutation frequencies by state
plot_mut_num <- ggplot(mut_counts_freqs, aes(fill=TREATMENT, y=MU_COUNT, x=TREATMENT)) +
  geom_boxplot() +
  xlab("") + 
  ylab("Mutation Counts") +
  ggtitle("Mutation Counts") +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, vjust=0.9, position = position_dodge(width = 1), size = 3 ) +
  facet_grid(cols=vars(POPULATION), scales = "free", drop = T) +
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1),legend.position = "none")
ggsave(plot=plot_mut_num, filename = paste0(mutation_dir,"/Mutation_count_state.svg"), device = "svg", width = 25, height = 10, units = "cm")
ggsave(plot=plot_mut_num, filename = paste0(mutation_dir,"/Mutation_count_state.png"), device = "png", width = 25, height = 10, units = "cm")

