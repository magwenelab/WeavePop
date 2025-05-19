log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(patchwork))

print("Reading files...")
metadata <- read.csv(snakemake@input[[1]], header = TRUE, stringsAsFactors = TRUE)
chrom_names <- read.csv(snakemake@input[[2]], header = TRUE, colClasses = "factor")
map_stats <- read.table(snakemake@input[[3]], header = TRUE, stringsAsFactors = TRUE, sep = "\t")

print("Joining and arranging data...")
metadata <- metadata %>%
    select(sample, strain, lineage) %>%
    mutate(name = paste(strain, sample, sep = " "))

chrom_names <- select(chrom_names, lineage, accession, chromosome)

map_stats <- left_join(metadata, map_stats, by = "sample")
map_stats$name <- reorder(map_stats$name, -map_stats$genome_mean_depth_good, sum)

depth_good <- map_stats %>%
    select(sample, name, lineage, Mean = genome_mean_depth_good, Median = genome_median_depth_good) %>%
    pivot_longer(cols = c(Mean, Median), names_to = "measurement", values_to = "value") %>%
    mutate(quality = "Good quality mappings")
depth_raw <- map_stats %>%
    select(sample, name, lineage, Mean = genome_mean_depth_raw, Median = genome_median_depth_raw) %>%
    pivot_longer(cols = c(Mean, Median), names_to = "measurement", values_to = "value") %>%
    mutate(quality = "All mappings")
depth <- rbind(depth_good, depth_raw)

print("Getting plot parameters...")
topylim <- max(depth$value) + max(depth$value/ 10)
raw_color = "gray50"
good_color = "black" 
color_quality = c("Good quality mappings" = good_color, "All mappings" = raw_color)
shape_stat <- c("Mean" = 16, "Median" = 15)


print("Plotting genome-wide read depth...")
g <- ggplot(depth) +
    geom_point(aes(x = name, y = value, shape = measurement, color = quality)) +
    ylim(0, topylim) +
    facet_grid(~lineage, scale = "free_x", space = "free_x") +
    scale_color_manual(name= "", values= color_quality)+
    scale_shape_manual(values = c(16,15, 17), name = NULL)+
    theme_bw() +
    theme(panel.background = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())+
    labs(title = "Genome-Wide Read Depth",
         y = "Read Depth (X)",
         x = "")

print("Joining and arranging data...")

stats_metad <- map_stats %>%
    select(sample, name, lineage, strain, percent_only_mapped, percent_unmapped, percent_properly_paired, percent_low_mapq, percent_inter_mapq, percent_high_mapq)

stats_long <- stats_metad %>%
    pivot_longer(cols = -c(sample, name, lineage, strain), names_to = "metric", values_to = "value")

stats_reads <- stats_long %>%
    filter(metric %in% c("percent_only_mapped", "percent_unmapped", "percent_properly_paired"))
stats_reads$metric <- factor(stats_reads$metric, levels = c("percent_unmapped", "percent_only_mapped", "percent_properly_paired"),
                              labels = c("Unmapped", "Mapped", "Mapped and properly paired"))

stats_qualit <- stats_long %>%
    filter(metric %in% c("percent_low_mapq", "percent_inter_mapq", "percent_high_mapq"))
stats_qualit$metric <- factor(stats_qualit$metric, levels = c("percent_low_mapq","percent_inter_mapq","percent_high_mapq"),
                              labels = c("Low MAPQ", "Intermediate MAPQ", "High MAPQ"))

print("Getting plot parameters...")
palette_reads <- brewer.pal(n = length(unique(stats_reads$metric)), name = "BuPu")
palette_qualit <- brewer.pal(n = length(unique(stats_qualit$metric)), name = "BuGn")

print("Plotting percentage of reads by mapping status...")
reads <- ggplot()+
    geom_bar(data = stats_reads, aes(x = name, y = value, fill = metric), stat = "identity")+
    facet_grid(~ lineage, scales = "free", space = "free_x")+
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    labs(x = "", y = "Percentage of Reads", fill = "Metric", title = "Percentage of Reads by Mapping Status")+
    scale_fill_manual(values = palette_reads, name = "")

print("Plotting percentage of reads by mapping quality...")
mapq <- ggplot() +
    geom_bar(data = stats_qualit, aes(x = name, y = value, fill = metric), stat = "identity") +
    facet_grid(~ lineage, scales = "free", space = "free_x") +
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5)) +
    labs(x = "", y = "Percentage of Reads", fill = "Metric", title = "Percentage of Mapped Reads by Mapping Quality") +
    scale_fill_manual(values = palette_qualit, name = "")

print("Joining plots...")
plot <- g/reads/mapq

print("Saving plot...")
gscale = snakemake@params[[1]]
ggsave(snakemake@output[[1]], plot = plot, units = "in", height = 9, width = 16, scale = gscale)
print("Done!")