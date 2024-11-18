log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(RColorBrewer))

print("Reading files...")
raw_stats_chroms <- read.delim(snakemake@input[[1]], sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
good_stats_chroms <- read.delim(snakemake@input[[2]], sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
chrom_names <- read.csv(snakemake@input[[3]], sep = ",", header = FALSE, col.names = c("lineage", "accession", "chromosome"), stringsAsFactors = TRUE, na = c("", "N/A"))
sample <- unique(good_stats_chroms$sample)

print("Filtering chromosome names...")
chrom_names <- chrom_names %>%
  filter(accession %in% unique(good_stats_chroms$accession) )

print("Ordering chromosome names...")
chrom_names['accession_chromosome'] <- paste(chrom_names$chromosome, chrom_names$accession, sep = "xxx")
unique_levels <- unique(chrom_names$accession_chromosome)
chrom_names$accession_chromosome <- factor(chrom_names$accession_chromosome, levels = unique_levels)

print("Joining and pivoting data...")
good_stats_chroms <- left_join(good_stats_chroms, chrom_names, by = "accession")
raw_stats_chroms <- left_join(raw_stats_chroms, chrom_names, by = "accession")

good_stats_long <- good_stats_chroms %>%
  pivot_longer(c(chrom_mean, chrom_median), names_to = "measurement", values_to = "value")
raw_stats_long <- raw_stats_chroms %>%
  pivot_longer(c(chrom_mean, chrom_median), names_to = "measurement", values_to = "value")

print("Getting plot parameters...")
toplim <- max(raw_stats_long$value) + max(raw_stats_long$value)/10
lineage <- unique(good_stats_chroms$lineage)
raw_color = "gray50"
good_color = "black" 
color_quality = c("Good quality alignments" = good_color, "All alignments" = raw_color)

print("Ploting good quality Read depth...")            
plot <- ggplot()+
  ylim(0,toplim) +
  geom_hline(aes(yintercept = unique(raw_stats_long$global_median),linetype = "global median", color= "All alignments"))+
  geom_hline(aes(yintercept = unique(raw_stats_long$global_mean), linetype = "global mean", color= "All alignments"))+
  geom_hline(aes(yintercept = unique(good_stats_long$global_median),linetype = "global median", color = "Good quality alignments"))+
  geom_hline(aes(yintercept = unique(good_stats_long$global_mean), linetype = "global mean", color = "Good quality alignments"))+
  geom_hline(aes(yintercept = unique(good_stats_long$global_mode), linetype = "global mode", color = "Good quality alignments"))+
  geom_point(data = raw_stats_long, aes(x = accession_chromosome, y = value, shape = measurement, color= "All alignments"))+ 
  geom_point(data = good_stats_long, aes(x = accession_chromosome, y = value, shape = measurement, color = "Good quality alignments"))+ 
  labs(y = "Read depth", x = "Chromosome", title = paste("Lineage:", lineage," Sample:", sample,  sep = " "))+
  theme_bw()+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_shape_manual(values = c(16,15), name = NULL, labels = c("Mean", "Median"))+
  scale_linetype_manual(values = c("dashed","dotted", "solid"), name = NULL, labels =c("Global mean", "Global median", "Global mode"))+
  scale_color_manual(name= "Alignment quality", values= color_quality)+
  scale_x_discrete(labels = function(x) gsub("xxx.*", "", x))

print("Saving plot...")
ggsave(snakemake@output[[1]], plot = plot, units = "in", height = 7.5, width = 7, dpi = 600 )
print("Done!") 