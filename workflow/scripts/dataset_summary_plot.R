log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggnewscale))

print("Reading input parameters...")
metadata_input <- snakemake@input$metadata
chrom_names_input <- snakemake@input$chroms
map_stats_input <- snakemake@input$stats

output <- snakemake@output$plot


# metadata_input <- "/PhastData/czirion/WeavePop_workdir/Cdeneoformans_global/NC/results/02.Dataset/metadata.csv"
# chrom_names_input <- "/PhastData/czirion/WeavePop_workdir/Cdeneoformans_global/NC/results/02.Dataset/chromosomes.csv"
# map_stats_input <- "/PhastData/czirion/WeavePop_workdir/Cdeneoformans_global/NC/results/02.Dataset/depth_quality/mapping_stats.tsv"

# output <- "plot.png"

print("Reading files...")
metadata <- read.csv(metadata_input, header = TRUE, stringsAsFactors = TRUE)
chrom_names <- read.csv(chrom_names_input, header = TRUE, colClasses = "factor")
map_stats <- read.table(map_stats_input, header = TRUE, stringsAsFactors = TRUE, sep = "\t")

# Add column to show red ticks in variables with warning
map_stats <- map_stats %>%
    mutate(depth_ticks = ifelse(grepl('Depth',quality_warning), "w", NA),
            mapq_ticks = ifelse(grepl('MAPQ',quality_warning), "w", NA),
            mapped_ticks = ifelse(grepl('Mapped',quality_warning), "w", NA),
            cov_ticks = ifelse(grepl('Coverage',quality_warning), "w", NA))

map_stats$quality_warning <- ifelse(map_stats$quality_warning == "", NA, map_stats$quality_warning)

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

map_stats <- left_join(metadata, map_stats, by = "sample")%>%
    arrange(name)%>%
    mutate(name = factor(name, levels = unique(name)),
            color = ifelse(is.na(quality_warning), "black", "red"),
            colored_label = paste0("<span style='color:", color, "'>", name, "</span>"))%>%
    mutate(colored_label = factor(colored_label, levels = unique(colored_label)))
# map_stats$name <- reorder(map_stats$name, -map_stats$genome_mean_depth_good, sum)

depth_good <- map_stats %>%
    select(sample, name, ref_genome, depth_ticks, Mean = genome_mean_depth_good, Median = genome_median_depth_good) %>%
    pivot_longer(cols = c(Mean, Median), names_to = "measurement", values_to = "value") %>%
    mutate(quality = "Good quality reads")
depth_raw <- map_stats %>%
    select(sample, name, ref_genome, depth_ticks, Mean = genome_mean_depth_raw, Median = genome_median_depth_raw) %>%
    pivot_longer(cols = c(Mean, Median), names_to = "measurement", values_to = "value") %>%
    mutate(quality = "All reads")
depth <- rbind(depth_good, depth_raw)

coverage_good <- map_stats %>%
    select(sample, name, ref_genome,cov_ticks, Coverage = coverage_good) %>%
    pivot_longer(cols = Coverage, names_to = "measurement", values_to = "value") %>%
    mutate(quality = "Good quality reads")
coverage_raw <- map_stats %>%
    select(sample, name, ref_genome, cov_ticks, Coverage = coverage_raw) %>%
    pivot_longer(cols = Coverage, names_to = "measurement", values_to = "value") %>%
    mutate(quality = "All reads")
coverage <- rbind(coverage_good, coverage_raw)


print("Getting plot parameters...")
topylim <- max(depth$value) + max(depth$value/ 10)
raw_color = "gray50"
good_color = "black" 
color_quality = c("Good quality reads" = good_color, "All reads" = raw_color)
shape_stat <- c("Mean" = 16, "Median" = 15)

print("Plotting genome-wide read depth...")
g <- ggplot(depth) +
    geom_point(aes(x = name, y = value, shape = measurement, color = quality)) +
    ylim(0, topylim) +
    facet_grid(~ref_genome, scale = "free_x", space = "free_x") +
    scale_color_manual(name= "", values= color_quality)+
    scale_shape_manual(values = c(16,15,17), name = NULL)+
    new_scale_color()+
    new_scale("shape") +
    geom_point(aes(x = name, y = 0, shape = depth_ticks), color = "red") +
    scale_shape_manual(values = c(8), name = NULL)+
    guides(shape = "none") +
    theme_bw() +
    theme(panel.background = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
          axis.text.x = element_blank(),
          #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5), 
          axis.ticks.x = element_blank())+
    labs(title = "Read Depth",
         y = "Depth (X)",
         x = "")
print("Plotting coverage...")
min_cov <- as.integer(round(min(coverage$value, na.rm = TRUE),0) - 2)
c <- ggplot(coverage) +
    geom_point(aes(x = name, y = value, color = quality)) +
    facet_grid(~ref_genome, scale = "free_x", space = "free_x") +
    scale_color_manual(name= "", values= color_quality)+
    scale_shape_manual(values = c(16,15, 17), name = NULL)+
    new_scale_color()+
    new_scale("shape") +
    geom_point(aes(x = name, y = min_cov, shape = cov_ticks), color = "red") +
    scale_shape_manual(values = c(8), name = NULL)+
    scale_y_continuous(breaks = pretty_breaks())+
    guides(shape = "none") +
    theme_bw() +
    theme(panel.background = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
          axis.text.x = element_blank(),
          #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5), 
          axis.ticks.x = element_blank())+
    labs(title = "Coverage",
         y = "Coverage",
         x = "")

print("Joining and arranging data...")

stats_metad <- map_stats %>%
    select(sample, name,colored_label, ref_genome, quality_warning,
        mapped_ticks, mapq_ticks,
        percent_only_mapped, percent_unmapped,percent_properly_paired,
        percent_low_mapq, percent_inter_mapq, percent_high_mapq,
        percent_filtered_vars, percent_removed_vars,
        n_raw_vars)

stats_long <- stats_metad %>%
    pivot_longer(cols = -c(sample, name, colored_label, ref_genome, quality_warning, mapped_ticks, mapq_ticks), names_to = "metric", values_to = "value")

stats_reads <- stats_long %>%
    filter(metric %in% c("percent_only_mapped", "percent_unmapped", "percent_properly_paired"))
stats_reads$metric <- factor(stats_reads$metric, levels = c("percent_unmapped", "percent_only_mapped", "percent_properly_paired"),
                              labels = c("Unmapped", "Mapped", "Mapped and\nproperly paired"))

stats_qualit <- stats_long %>%
    filter(metric %in% c("percent_low_mapq", "percent_inter_mapq", "percent_high_mapq"))
stats_qualit$metric <- factor(stats_qualit$metric, levels = c("percent_low_mapq","percent_inter_mapq","percent_high_mapq"),
                              labels = c("Low MAPQ", "Intermediate MAPQ", "High MAPQ"))

stats_vars <- stats_long %>%
    filter(metric %in% c("percent_filtered_vars", "percent_removed_vars"))
stats_vars$metric <- factor(stats_vars$metric, levels = c("percent_removed_vars","percent_filtered_vars"),
                              labels = c("Removed variants","Filtered variants"))

stats_total_vars <- stats_long %>%
    filter(metric %in% c("n_raw_vars"))
stats_total_vars$metric <- factor(stats_total_vars$metric, levels = c("n_raw_vars"),
                              labels = c("Total variants"))

print("Getting plot parameters...")
palette_reads <- brewer.pal(n = length(unique(stats_reads$metric)), name = "BuPu")
palette_qualit <- brewer.pal(n = length(unique(stats_qualit$metric)), name = "BuGn")
palette_vars <- brewer.pal(n = 11, name = "BrBG")[c(5,3)]

print("Plotting percentage of reads by mapping status...")
reads <- ggplot(data = stats_reads)+
    geom_bar(aes(x = name, y = value, fill = metric), stat = "identity")+
    facet_grid(~ ref_genome, scales = "free", space = "free_x")+
    new_scale_color()+
    geom_point(aes(x = name, y = 0, shape = mapped_ticks), color = "red") +
    scale_shape_manual(values = c(8), name = NULL)+
    guides(shape = "none") +
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
          axis.text.x = element_blank(),
          #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5), 
          axis.ticks.x = element_blank())+
    labs(x = "", y = "% Reads", fill = "Metric", title = "Percentage of Reads by Mapping Status")+
    scale_fill_manual(values = palette_reads, name = "")

print("Plotting percentage of reads by mapping quality...")
mapq <- ggplot(data = stats_qualit) +
    geom_bar(aes(x = name, y = value, fill = metric), stat = "identity") +
    facet_grid(~ ref_genome, scales = "free", space = "free_x") +
    new_scale_color()+
    geom_point(aes(x = name, y = 0, shape = mapq_ticks), color = "red") +
    scale_shape_manual(values = c(8), name = NULL)+
    guides(shape = "none") +
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
          #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5), 
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    labs(x = "", y = "% Reads", fill = "Metric", title = "Percentage of Mapped Reads by Mapping Quality") +
    scale_fill_manual(values = palette_qualit, name = "")

print("Plotting absolute of removed and filtered variants...")
vars_total <- ggplot() +
    geom_point(data = stats_total_vars, aes(x = name, y = value, fill = metric), stat = "identity") +
    facet_grid(~ ref_genome, scales = "free", space = "free_x") +
    theme_bw() +
    theme(panel.background = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
          #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5), 
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    labs(x = "", y = "# Variants", fill = "", title = "Number of Variants")+
    scale_y_continuous(labels = scales::comma)

print("Plotting percentage of removed and filtered variants...")
vars <- ggplot() +
    geom_bar(data = stats_vars, aes(x = colored_label, y = value, fill = metric), stat = "identity") +
    facet_grid(~ ref_genome, scales = "free", space = "free_x") +
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
          axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1, vjust = 0.5, size = 5)) +
    labs(x = "", y = "% Variants", fill = "Metric", title = "Percentage of Variants by Filtered Status") +
    scale_fill_manual(values = palette_vars, name = "")



print("Joining plots...")
plot <- g/c/reads/mapq/vars_total/vars

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

gheight <- 13

ggsave(output, plot = plot, units = "in", height = gheight, width = gwidth, scale = gscale, limitsize = FALSE)

print("Done!")