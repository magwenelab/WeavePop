log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(RColorBrewer))

print("Loading and processing data")

# raw_stats_chroms <- read.delim("/FastData/czirion/ashton/results/samples/mosdepth/ERS542301/depth_by_chrom_raw.tsv", sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
raw_stats_chroms <- read.delim(snakemake@input[[1]], sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))

# good_stats_chroms <- read.delim("/FastData/czirion/ashton/results/samples/mosdepth/ERS542301/depth_by_chrom_good.tsv", sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
good_stats_chroms <- read.delim(snakemake@input[[2]], sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))

# chrom_names <- read.csv("/FastData/czirion/ashton/config/chromosome_names.csv", sep = ",", header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"), stringsAsFactors = TRUE, na = c("", "N/A"))
chrom_names <- read.csv(snakemake@input[[3]], sep = ",", header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"), stringsAsFactors = TRUE, na = c("", "N/A"))
sample <- unique(good_stats_chroms$Sample)

chrom_names <- chrom_names %>%
  filter(Accession %in% unique(good_stats_chroms$Accession) )

chrom_names['Accession_Chromosome'] <- paste(chrom_names$Chromosome, chrom_names$Accession, sep = "xxx")
unique_levels <- unique(chrom_names$Accession_Chromosome)
chrom_names$Accession_Chromosome <- factor(chrom_names$Accession_Chromosome, levels = unique_levels)


print("Joining and pivoting data")
good_stats_chroms <- left_join(good_stats_chroms, chrom_names, by = "Accession")
raw_stats_chroms <- left_join(raw_stats_chroms, chrom_names, by = "Accession")

good_stats_long <- good_stats_chroms %>%
  pivot_longer(c(Chrom_Mean, Chrom_Median), names_to = "Measurement", values_to = "Value")
raw_stats_long <- raw_stats_chroms %>%
  pivot_longer(c(Chrom_Mean, Chrom_Median), names_to = "Measurement", values_to = "Value")

print("Making plot parameters")
toplim <- max(raw_stats_long$Value) + max(raw_stats_long$Value)/10
lineage <- unique(good_stats_chroms$Lineage)

raw_color = "#B3B3B3"
good_color = "#666666" 
color_quality = c("Good quality alignments" = good_color, "All alignments" = raw_color)

print("Ploting good quality Read depth")            

plot <- ggplot()+
  ylim(0,toplim) +
  geom_hline(aes(yintercept = unique(raw_stats_long$Global_Median),linetype = "Global median", color= "All alignments"))+
  geom_hline(aes(yintercept = unique(raw_stats_long$Global_Mean), linetype = "Global mean", color= "All alignments"))+
  geom_hline(aes(yintercept = unique(good_stats_long$Global_Median),linetype = "Global median", color = "Good quality alignments"))+
  geom_hline(aes(yintercept = unique(good_stats_long$Global_Mean), linetype = "Global mean", color = "Good quality alignments"))+
  geom_hline(aes(yintercept = unique(good_stats_long$Global_Mode), linetype = "Global mode", color = "Good quality alignments"))+
  geom_point(data = raw_stats_long, aes(x = Accession_Chromosome, y = Value, shape = Measurement, color= "All alignments"))+ 
  geom_point(data = good_stats_long, aes(x = Accession_Chromosome, y = Value, shape = Measurement, color = "Good quality alignments"))+ 
  labs(y = "Read depth", x = "Chromosome", title = paste("Lineage:", lineage," Sample:", sample,  sep = " "))+
  theme_bw()+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_shape_manual(values = c(16,15), name = NULL, labels = c("Mean", "Median"))+
  scale_linetype_manual(values = c("dashed","dotted", "solid"), name = NULL)+
  scale_color_manual(name= "Alignment quality", values= color_quality)

plot <- plot +
  scale_x_discrete(labels = function(x) gsub("xxx.*", "", x))

print("Saving plot")
ggsave(snakemake@output[[1]], plot = plot, units = "in", height = 7.5, width = 7, dpi = 600 )
print("Done!") 