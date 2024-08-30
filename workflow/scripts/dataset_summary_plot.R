log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(patchwork))

print("Reading files")
#metadata <- read.csv("config/sample_metadata_filtered.csv", header = TRUE, stringsAsFactors = TRUE)
# chrom_names <- read.csv("config/chromosome_names.csv", header = TRUE, col.names = c("lineage", "Accession", "Chromosome"), colClasses = "factor")
# good_stats <- read.delim("results/dataset/files/depth_by_chrom_good.tsv", sep = "\t", header = TRUE, stringsAsFactors = TRUE)
# raw_stats <- read.delim("results/dataset/files/depth_by_chrom_raw.tsv", sep = "\t", header = TRUE, stringsAsFactors = TRUE)

metadata <- read.csv(snakemake@input[[1]], header = TRUE, stringsAsFactors = TRUE)
metadata <- metadata %>%
    select(sample, strain, lineage) %>%
    mutate(name = paste(strain, sample, sep = " "))

chrom_names <- read.csv(snakemake@input[[2]], header = TRUE, col.names = c("lineage", "Accession", "Chromosome"), colClasses = "factor")
chrom_names <- select(chrom_names, Accession, Chromosome)

good_stats <- read.delim(snakemake@input[[3]], sep = "\t", header = TRUE, stringsAsFactors = TRUE)
good_stats <- rename(good_stats, sample = Sample)
good_stats <- left_join(good_stats, metadata, by = "sample")
good_stats <- left_join(good_stats, chrom_names, by = "Accession")

raw_stats <- read.delim(snakemake@input[[4]], sep = "\t", header = TRUE, stringsAsFactors = TRUE)
raw_stats <- rename(raw_stats, sample = Sample)
raw_stats <- left_join(raw_stats, metadata, by = "sample")
raw_stats <- left_join(raw_stats, chrom_names, by = "Accession")

print("Making plot parameters")
topylim <- max(good_stats$Global_Mean) + max(good_stats$Global_Mean / 10)
raw_color = "#B3B3B3"
good_color = "#666666" 
color_quality = c("Good quality mappings" = good_color, "All mappings" = raw_color)
shape_stat <- c("Mean" = 16, "Median" = 15)

print("Binding and pivoting data")
all <- bind_rows(good_stats %>% mutate(Quality = "Good quality mappings"), raw_stats %>% mutate(Quality = "All mappings"))
all$name <- reorder(all$name, -all$Global_Mean, sum)
all <- all %>%
        select(sample, name, Mean = Global_Mean, Median = Global_Median, Mode = Global_Mode, Quality, lineage)%>%
        pivot_longer(cols = c(Mean, Median, Mode), names_to = "Measurement", values_to = "Value")%>%
        as.data.frame()

print("Plotting genome-wide read depth")
topylim <- max(all$Value) + max(all$Value) / 10
g <- ggplot(all) +
    geom_point(aes(x = name, y = Value, shape = Measurement, color = Quality)) +
    ylim(0, topylim) +
    facet_grid(~lineage, scale = "free_x", space = "free_x") +
    scale_color_manual(name= "", values= color_quality)+
    scale_shape_manual(values = c(16,15, 17), name = NULL, labels = c("Mean", "Median", "Mode"))+
    theme_bw() +
    theme(panel.background = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())+
    labs(title = "Genome-wide read depth per sample",
         y = "Read depth (X)",
         x = "")

print("Reading and processing files")
# map_stats <- read.table("results/dataset/files/mapping_stats.tsv",header = TRUE, stringsAsFactors = TRUE, sep = "\t")
map_stats <- read.table(snakemake@input[[5]], header = TRUE, stringsAsFactors = TRUE, sep = "\t")

map_stats_metad <- merge(map_stats, metadata, by = "sample")
stats_metad <- map_stats_metad %>%
    select(sample, name, lineage, strain, percent_only_mapped, percent_unmapped, percent_properly_paired, percent_low_mapq, percent_inter_mapq, percent_high_mapq)

stats_long <- stats_metad %>%
    pivot_longer(cols = -c(sample, name, lineage, strain), names_to = "Metric", values_to = "Value")
stats_long$name <- factor(stats_long$name, levels = levels(all$name))

stats_reads <- stats_long %>%
    filter(Metric %in% c("percent_only_mapped", "percent_unmapped", "percent_properly_paired"))
stats_reads$Metric <- factor(stats_reads$Metric, levels = c("percent_unmapped", "percent_only_mapped", "percent_properly_paired"),
                              labels = c("Unmapped", "Mapped", "Mapped and properly paired"))

stats_qualit <- stats_long %>%
    filter(Metric %in% c("percent_low_mapq", "percent_inter_mapq", "percent_high_mapq"))
stats_qualit$Metric <- factor(stats_qualit$Metric, levels = c("percent_low_mapq","percent_inter_mapq","percent_high_mapq"),
                              labels = c("Low MAPQ", "Intermediate MAPQ", "High MAPQ"))

print("Making plot parameters")
palette_reads <- brewer.pal(n = length(unique(stats_reads$Metric)), name = "BuPu")
palette_qualit <- brewer.pal(n = length(unique(stats_qualit$Metric)), name = "BuGn")

print("Plotting percentage of reads by mapping status")
reads <- ggplot()+
    geom_bar(data = stats_reads, aes(x = name, y = Value, fill = Metric), stat = "identity")+
    facet_grid(~ lineage, scales = "free", space = "free_x")+
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    labs(x = "", y = "Percentage of reads", fill = "Metric", title = "Percentage of reads by mapping status")+
    scale_fill_manual(values = palette_reads, name = "")

print("Plotting percentage of reads by mapping quality")
mapq <- ggplot() +
    geom_bar(data = stats_qualit, aes(x = name, y = Value, fill = Metric), stat = "identity") +
    facet_grid(~ lineage, scales = "free", space = "free_x") +
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5)) +
    labs(x = "", y = "Percentage of reads", fill = "Metric", title = "Percentage of mapped reads by mapping quality") +
    scale_fill_manual(values = palette_qualit, name = "")

print("Joining plots")
plot <- g/reads/mapq

print("Saving plot")
gscale = snakemake@params[[1]]
ggsave(snakemake@output[[1]], plot = plot, units = "in", height = 9, width = 16, scale = gscale)
print("Done!")