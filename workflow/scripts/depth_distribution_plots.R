log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scales))

print("Reading files and joining data with chromosome names")
# setwd("./results/samples/samtools/SRS8318899")
# cov<- read.delim("distrib_cov.tsv", header = TRUE, stringsAsFactors = TRUE)
# chrom_names <- read.csv("../../../../config/chromosome_names.csv", header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))
# sample <- "/mnt/FastData/czirion/DiversityPipeline/results/samples/samtools/SRS8318899/distrib_cov.csv"
sample <- snakemake@wildcards$sample
cov<- read.table(snakemake@input[[1]], header = TRUE, stringsAsFactors = TRUE, sep = "\t")
global_mode <- read.delim(snakemake@input[[2]], header = TRUE, stringsAsFactors = TRUE)
print(global_mode)
chrom_names <- read.csv(snakemake@input[[3]], header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))

cov[is.na(cov)] <- 0
cov <- left_join(cov, chrom_names, by = "Accession")


print("Plotting depth distribution by chromosome")
raw_color = "#B3B3B3"
good_color = "#666666" 
color_quality = c("Good quality alignments" = good_color, "All alignments" = raw_color)
lineage <- levels(as.factor(cov$Lineage))

plot_chrom <- ggplot(cov, aes(x=Depth))+
  geom_col(aes(y = Count_Raw, fill= "All alignments"))+ 
  geom_col(aes(y = Count_Good, fill= "Good quality alignments"))+ 
  facet_wrap(~Chromosome,ncol = 2)+
  scale_y_log10(name = "Number of Sites", labels = comma)+
  scale_x_continuous(name = "Depth (X) ", labels = comma, n.breaks = 10)+
  scale_fill_manual(name= "Alignment quality", values= color_quality)+
  theme(legend.position="none")+
  labs(title = paste("Lineage:",lineage, " Sample:", sample,  sep = " "))+
  theme_bw()

pheight <- 0.5 + length(unique((cov$Chromosome)))/2
pwidth <- pheight * 1.78

cov_global <- cov %>%
  select(Depth, Count_Good, Count_Raw)%>%
  group_by(Depth)%>%
  summarize(Count_Good_Global = sum(Count_Good), Count_Raw_Global = sum(Count_Raw))%>%
  ungroup()

print("Plotting genome-wide depth distribution")
plot_global <- ggplot()+
  geom_col(data = cov_global, aes(x=Depth, y = Count_Raw_Global, fill= "All alignments"))+ 
  geom_col(data = cov_global, aes(x=Depth,y = Count_Good_Global, fill= "Good quality alignments"))+ 
  geom_vline(data = global_mode, aes(xintercept = Global_Mode), color ="red")+
  scale_y_log10(name = "Number of Sites", labels = comma)+
  scale_x_continuous(name = "Depth (X) ", labels = comma, n.breaks = 10)+
  scale_fill_manual(name= "Alignment quality", values= color_quality)+
  theme(legend.position="none")+
  labs(title = paste("Lineage:",lineage, " Sample:", sample,  sep = " "))+
  theme_bw()

print("Saving plots")
ggsave(snakemake@output[[1]], plot = plot_chrom, units = "in", height = pheight, width = pwidth, dpi = 600)
ggsave(snakemake@output[[2]], plot = plot_global,  units = "in", height = 4.5, width = 8, dpi = 600)
print("Done!")
