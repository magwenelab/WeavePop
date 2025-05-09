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
chrom_names <- read.csv(snakemake@input[[2]], header = TRUE, stringsAsFactors = TRUE)
chrom_depth <- read.delim(snakemake@input[[3]], sep = "\t", header = TRUE, stringsAsFactors = TRUE)

print("Filtering chromosome names...")
chrom_names <- chrom_names %>%
    select(lineage, accession, chromosome)

print("Joining and arranging data...")
metadata <- mutate(metadata, name = paste(strain, sample, sep = " "))
chrom_depth <- left_join(chrom_depth, metadata, by = "sample")
chrom_depth <- left_join(chrom_depth, chrom_names, by = c("accession", "lineage"))

print("Getting plot parameters...")
toplim <- ceiling(max(chrom_depth$norm_chrom_median))
values <- seq(0, toplim, by = 1)
gscale = snakemake@params[[2]]

print("Plotting...")
medianplot <- ggplot(chrom_depth, aes(x = reorder(name, -genome_median_depth, sum), y = norm_chrom_median)) +
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
    labs(title = "Normalized Median of the Depth of the Windows Along Each Chromosome",
         x = "Sample",
         y = "Normalized Depth")

print("Saving plot...")
ggsave(snakemake@output[[1]], plot = medianplot, units = "in", height = 9, width = 16, scale = gscale)
print("Done!")