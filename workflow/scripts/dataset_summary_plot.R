log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(patchwork))

print("Reading input parameters...")
metadata_input <- snakemake@input$metadata
chrom_names_input <- snakemake@input$chroms
map_stats_input <- snakemake@input$stats

output <- snakemake@output$plot

print("Reading files...")
metadata <- read.csv(metadata_input, header = TRUE, stringsAsFactors = TRUE)
chrom_names <- read.csv(chrom_names_input, header = TRUE, colClasses = "factor")
map_stats <- read.table(map_stats_input, header = TRUE, stringsAsFactors = TRUE, sep = "\t")

print("Joining and arranging data...")

if ("strain" %in% colnames(metadata)){
    metadata <- metadata %>%
        select(sample, strain, ref_genome) %>%
        mutate(name = paste(strain, sample, sep = " ")) 
    } else {
    metadata <- metadata %>%
        select(sample, ref_genome) %>%
        mutate(name = sample)}

chrom_names <- select(chrom_names, ref_genome, accession, chromosome)

map_stats <- left_join(metadata, map_stats, by = "sample")
map_stats$name <- reorder(map_stats$name, -map_stats$genome_mean_depth_good, sum)

depth_good <- map_stats %>%
    select(sample, name, ref_genome, Mean = genome_mean_depth_good, Median = genome_median_depth_good) %>%
    pivot_longer(cols = c(Mean, Median), names_to = "measurement", values_to = "value") %>%
    mutate(quality = "Good quality mappings")
depth_raw <- map_stats %>%
    select(sample, name, ref_genome, Mean = genome_mean_depth_raw, Median = genome_median_depth_raw) %>%
    pivot_longer(cols = c(Mean, Median), names_to = "measurement", values_to = "value") %>%
    mutate(quality = "All mappings")
depth <- rbind(depth_good, depth_raw)

coverage_good <- map_stats %>%
    select(sample, name, ref_genome, Coverage = coverage_good) %>%
    pivot_longer(cols = Coverage, names_to = "measurement", values_to = "value") %>%
    mutate(quality = "Good quality mappings")
coverage_raw <- map_stats %>%
    select(sample, name, ref_genome, Coverage = coverage_raw) %>%
    pivot_longer(cols = Coverage, names_to = "measurement", values_to = "value") %>%
    mutate(quality = "All mappings")
coverage <- rbind(coverage_good, coverage_raw)


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
    facet_grid(~ref_genome, scale = "free_x", space = "free_x") +
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

print("Plotting coverage...")
c <- ggplot(coverage) +
    geom_point(aes(x = name, y = value, color = quality)) +
    facet_grid(~ref_genome, scale = "free_x", space = "free_x") +
    scale_color_manual(name= "", values= color_quality)+
    scale_shape_manual(values = c(16,15, 17), name = NULL)+
    theme_bw() +
    theme(panel.background = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())+
    labs(title = "Coverage",
         y = "Percentage of Coverage",
         x = "")

print("Joining and arranging data...")

stats_metad <- map_stats %>%
    select(sample, name, ref_genome, strain, quality_warning,
        percent_only_mapped, percent_unmapped,percent_properly_paired,
        percent_low_mapq, percent_inter_mapq, percent_high_mapq)

stats_long <- stats_metad %>%
    pivot_longer(cols = -c(sample, name, ref_genome, strain, quality_warning), names_to = "metric", values_to = "value")

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
stats_qualit <- stats_qualit %>%
    mutate(color = ifelse(quality_warning == "", "black", "red"),
            colored_label = paste0("<span style='color:", color, "'>", name, "</span>"))%>%
    arrange(factor(name, levels = levels(stats_qualit$name)))
stats_qualit$colored_label <- factor(stats_qualit$colored_label, levels = unique(stats_qualit$colored_label))

print("Plotting percentage of reads by mapping status...")
reads <- ggplot()+
    geom_bar(data = stats_reads, aes(x = name, y = value, fill = metric), stat = "identity")+
    facet_grid(~ ref_genome, scales = "free", space = "free_x")+
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
    geom_bar(data = stats_qualit, aes(x = colored_label, y = value, fill = metric), stat = "identity") +
    facet_grid(~ ref_genome, scales = "free", space = "free_x") +
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
          axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1, vjust = 0.5, size = 5)) +
    labs(x = "", y = "Percentage of Reads", fill = "Metric", title = "Percentage of Mapped Reads by Mapping Quality") +
    scale_fill_manual(values = palette_qualit, name = "")

print("Joining plots...")
plot <- g/c/reads/mapq

print("Saving plot...")
n_samples <- nrow(metadata)

if (n_samples <= 200) {
    gscale <- 1
    gwidth <- 6 + n_samples * 0.03
} else if (n_samples > 200 & n_samples <= 400) {
    gwidth <- 16
    gscale <- 2
} else if ( n_samples > 400 & n_samples <= 1000) {
    gwidth <- 18
    gscale <- 2.5
} else {
    gwidth <- 20
    gscale <- 3
}

gheight <- 9

ggsave(output, plot = plot, units = "in", height = gheight, width = gwidth, scale = gscale, limitsize = FALSE)

print("Done!")