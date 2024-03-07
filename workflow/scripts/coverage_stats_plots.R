log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(RColorBrewer))

print("Loading and processing data")
good_stats_chroms <- read.delim(snakemake@input[[1]], sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
raw_stats_chroms <- read.delim(snakemake@input[[2]], sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
chrom_names <- read.csv(snakemake@input[[3]], sep = ",", header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"), stringsAsFactors = TRUE, na = c("", "N/A"))
lineage <- unique(chrom_names$Lineage)
sample <- unique(good_stats_chroms$Sample)

good_stats_chroms <- left_join(good_stats_chroms, chrom_names, by = "Accession")
raw_stats_chroms <- left_join(raw_stats_chroms, chrom_names, by = "Accession")
good_stats_long <- good_stats_chroms %>%
  pivot_longer(c(Chrom_Mean, Chrom_Median), names_to = "Measurement", values_to = "Value")
raw_stats_long <- raw_stats_chroms %>%
  pivot_longer(c(Chrom_Mean, Chrom_Median), names_to = "Measurement", values_to = "Value")
toplim <- max(raw_stats_long$Value) + max(raw_stats_long$Value)/10


print("Ploting good quality coverage")            
raw_color = "#B3B3B3"
good_color = "#666666" 
color_quality = c("Good quality alignments" = good_color, "All alignments" = raw_color)

plot <- ggplot()+
  ylim(0,toplim) +
  geom_hline(aes(yintercept = unique(raw_stats_long$Global_Median),linetype = "Global median", color= "All alignments"))+
  geom_hline(aes(yintercept = unique(raw_stats_long$Global_Mean), linetype = "Global mean", color= "All alignments"))+
  geom_hline(aes(yintercept = unique(good_stats_long$Global_Median),linetype = "Global median", color = "Good quality alignments"))+
  geom_hline(aes(yintercept = unique(good_stats_long$Global_Mean), linetype = "Global mean", color = "Good quality alignments"))+
  geom_point(data = raw_stats_long, aes(x = factor(Chromosome, levels = as.character(sort(unique(Chromosome)))), y = Value, shape = Measurement, color= "All alignments"))+ 
  geom_point(data = good_stats_long, aes(x = factor(Chromosome, levels = as.character(sort(unique(Chromosome)))), y = Value, shape = Measurement, color = "Good quality alignments"))+ 
  labs(y = "Coverage", x = "Chromosome", title = paste(lineage, sample,  sep = " "))+
  theme_bw()+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_shape_manual(values = c(16,15), name = NULL, labels = c("Mean", "Median"))+
  scale_linetype_manual(values = c("solid","dotted"), name = NULL)+
  scale_color_manual(name= "Alignment quality", values= color_quality)

print("Saving plot")
ggsave(snakemake@output[[1]], plot = plot, units = "in", height = 7.5, width = 7)
print("Done") 