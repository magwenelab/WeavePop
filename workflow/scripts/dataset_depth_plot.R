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

print("Reading files...")
metadata <- read.csv(metadata_input, header = TRUE, stringsAsFactors = TRUE)
chrom_names <- read.csv(chrom_names_input, header = TRUE, stringsAsFactors = TRUE)
chrom_depth <- read.delim(chrom_depth_input, sep = "\t", header = TRUE, stringsAsFactors = TRUE)

chrom_depth <- chrom_depth %>%
    select("sample", "accession", "genome_depth", "chrom_norm_depth") %>%
    distinct()

print("Filtering chromosome names...")
chrom_names <- chrom_names %>%
    select(ref_genome, accession, chromosome)

print("Joining and arranging data...")
metadata <- mutate(metadata, name = paste(strain, sample, sep = " "))
chrom_depth <- left_join(chrom_depth, metadata, by = "sample")
chrom_depth <- left_join(chrom_depth, chrom_names, by = c("accession", "ref_genome"))

print("Getting plot parameters...")
toplim <- ceiling(max(chrom_depth$chrom_norm_depth))
values <- seq(0, toplim, by = 1)

print("Plotting...")
medianplot <- ggplot(chrom_depth, aes(x = reorder(name, -genome_depth, sum), y = chrom_norm_depth)) +
    ylim(0, toplim) +
    facet_grid(scale = "free_x", space = "free_x", rows = vars(chromosome), cols = vars(ref_genome)) +
    theme_bw() +
    theme(panel.background = element_blank(), 
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5))+
    labs(title = "Normalized Median Depth per Chromosome",
         x = "",
         y = "Normalized Depth")

if (color_by %in% colnames(metadata)){
    medianplot <- medianplot +
        geom_col(aes(fill = get(color_by))) +
        scale_fill_brewer(palette = "Set2", name = str_to_title(color_by))
} else {
    medianplot <- medianplot +
        geom_col(fill = "black") 
}

n_samples <- length(unique(chrom_depth$sample))
n_chroms <- length(unique(chrom_depth$chromosome))

if (n_samples <= 200) {
    gscale <- 1
    gwidth <- 5 + n_samples * 0.03
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

gheight <- 4 + n_chroms * 0.5

print("Saving plot...")
ggsave(output, plot = medianplot, units = "in", height = gheight, width = gwidth, scale = gscale)
print("Done!")