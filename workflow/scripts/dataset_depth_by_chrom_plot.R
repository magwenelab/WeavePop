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
chrom_depth_input <- snakemake@input$cnv

output <- snakemake@output$plot

color_by <- snakemake@params$column
gscale <- snakemake@params$scale
gheight <- 9
gwidth <- 16


print("Reading files...")
metadata <- read.csv(metadata_input, header = TRUE, stringsAsFactors = TRUE)
chrom_names <- read.csv(chrom_names_input, header = TRUE, stringsAsFactors = TRUE)
chrom_depth <- read.delim(chrom_depth_input, sep = "\t", header = TRUE, stringsAsFactors = TRUE)

chrom_depth <- chrom_depth %>%
    select("sample", "accession", "genome_depth", "chrom_norm_depth") %>%
    distinct()

print("Filtering chromosome names...")
chrom_names <- chrom_names %>%
    select(lineage, accession, chromosome)

print("Joining and arranging data...")
metadata <- mutate(metadata, name = paste(strain, sample, sep = " "))
chrom_depth <- left_join(chrom_depth, metadata, by = "sample")
chrom_depth <- left_join(chrom_depth, chrom_names, by = c("accession", "lineage"))

print("Getting plot parameters...")
toplim <- ceiling(max(chrom_depth$chrom_norm_depth))
values <- seq(0, toplim, by = 1)

print("Plotting...")
medianplot <- ggplot(chrom_depth, aes(x = reorder(name, -genome_depth, sum), y = chrom_norm_depth)) +
    ylim(0, toplim) +
    facet_grid(scale = "free_x", space = "free_x", rows = vars(chromosome), cols = vars(lineage)) +
    theme_bw() +
    theme(panel.background = element_blank(), 
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5))+
    labs(title = "Normalized Median Depth per Chromosome",
         x = "Sample",
         y = "Normalized Depth")

if (color_by %in% colnames(metadata)){
    medianplot <- medianplot +
        geom_point(aes(color = get(color_by))) +
        scale_color_brewer(palette = "Set2", name = str_to_title(color_by))
} else {
    medianplot <- medianplot +
        geom_point(color = "black") 
}


print("Saving plot...")
ggsave(output, plot = medianplot, units = "in", height = gheight, width = gwidth, scale = gscale)
print("Done!")