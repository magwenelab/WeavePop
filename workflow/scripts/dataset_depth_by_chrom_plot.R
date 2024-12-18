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
chrom_names <- read.csv(snakemake@input[[2]], header = TRUE, col.names = c("lineage", "accession", "chromosome"), stringsAsFactors = TRUE)
good_stats <- read.delim(snakemake@input[[3]], sep = "\t", header = TRUE, stringsAsFactors = TRUE)

print("Filtering chromosome names...")
chrom_names <- chrom_names %>%
    select(accession, chromosome)

print("Joining and arranging data...")
metadata <- mutate(metadata, name = paste(strain, sample, sep = " "))
good_stats <- rename(good_stats, sample = sample)
good_stats <- left_join(good_stats, metadata, by = "sample")
good_stats <- left_join(good_stats, chrom_names, by = "accession")

print("Getting plot parameters...")
toplim <- ceiling(max(good_stats$norm_chrom_mean))
values <- seq(0, toplim, by = 1)
ylabel <- "Normalized depth"
gscale = snakemake@params[[2]]

print("Plotting...")
meanplot <- ggplot(good_stats, aes(x = reorder(name, -global_mean, sum), y = norm_chrom_mean)) +
    geom_point(aes(color = get(snakemake@params[[1]]))) +
    ylim(0, toplim) +
    facet_grid(scale = "free_x", space = "free_x", rows = vars(chromosome), cols = vars(lineage)) +
    scale_color_brewer(palette = "Set2", name = str_to_title(snakemake@params[[1]])) +
    theme_bw() +
    theme(panel.background = element_blank(), 
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5))+
    labs(title = "Normalized mean depth of chromosomes",
         x = "Sample",
         y = ylabel)

print("Saving plot...")
ggsave(snakemake@output[[1]], plot = meanplot, units = "in", height = 9, width = 16, scale = gscale)
print("Done!")