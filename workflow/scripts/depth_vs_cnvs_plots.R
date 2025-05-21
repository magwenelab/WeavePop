log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))

print("Reading files...")
cnv_chromosomes <- read.delim(snakemake@input[[1]], sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
metadata <- read.delim(snakemake@input[[2]], sep = ",", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
sample <- snakemake@wildcards$sample

print("Obtaining lineage of sample...")

lineage_name <- as.character(metadata$lineage[metadata$sample == sample])
strain_name <- as.character(metadata$strain[metadata$sample == sample])

chrom_metrics <- cnv_chromosomes %>%
    filter(cnv != "single_copy")

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
ggsave(snakemake@output[[1]], p,  width = 8, height = 6, dpi = 300)
print("Done.")