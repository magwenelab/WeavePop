log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))

print("Reading input parameters...")
cnv_chromosomes_input <- snakemake@input[[1]]
chrom_names_input <- snakemake@input[[2]]
metadata_input <- snakemake@input[[3]]

output <- snakemake@output[[1]]

sample <- snakemake@wildcards$sample

gscale <- 0.7
gheight <- 9
gwidth <- 16
gdpi <- 600

print("Reading files...")
cnv_chromosomes <- read.delim(cnv_chromosomes_input, sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
chrom_names <- read.csv(chrom_names_input,sep= ",", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
metadata <- read.delim(metadata_input, sep = ",", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))

print("Obtaining lineage of sample...")

lineage_name <- as.character(metadata$lineage[metadata$sample == sample])
strain_name <- as.character(metadata$strain[metadata$sample == sample])

chrom_metrics <- cnv_chromosomes %>%
    filter(cnv != "single_copy")

chrom_metrics <- chrom_metrics %>%
    left_join(chrom_names, by = "accession")

chrom_metrics$chromosome <- factor(chrom_metrics$chromosome)

   
print("Plotting...")
p <- ggplot(chrom_metrics, aes(x = coverage_percent, y = norm_chrom_median, color = chromosome, shape = cnv)) +
        geom_hline(yintercept = c(0, 1, 2), color = "black", linetype = "solid") +
        geom_point(size = 2) +
        geom_text_repel(aes(label = chromosome), size = 3, max.overlaps = 10, show.legend = FALSE) +
        scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
        theme_bw() +
        theme(legend.position = "right") +
        labs(title = "Normalized Depth vs.\nPercent of CNV Coverage per Chromosome",
            subtitle = paste("Lineage:", lineage_name, " Sample:", sample, "Strain:", strain_name, sep = " "),
             y = "Normalized Median Depth of Chromosome",
             x = "Percent of Chromosome Covered by CNVs",
             color = "Chromosome",
             shape = "CNV")

print("Saving plot...")
ggsave(output, p,  width = gwidth, height = gheight, dpi = gdpi, scale = gscale, units = "in")
print("Done.")