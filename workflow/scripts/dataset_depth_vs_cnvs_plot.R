log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(patchwork))

print("Reading files...")
depth_by_chrom <- read.delim(snakemake@input[[1]], sep = "\t", header = TRUE, stringsAsFactors = TRUE)
cnv_chroms <- read.delim(snakemake@input[[2]], sep = "\t", header = TRUE, stringsAsFactors = TRUE)
gscale = snakemake@params[[1]]


print("Joining and arranging data...")

chrom_metrics <- cnv_chroms %>%
    filter(cnv != "single_copy") %>%
    left_join(depth_by_chrom, by = c("sample", "accession"))%>%
    as.data.frame()

dups <- chrom_metrics %>%
    filter(cnv == "duplication")

dels <- chrom_metrics %>%
    filter(cnv == "deletion")

print("Plotting...")

dup <- ggplot(dups, aes(x = coverage_percent, y = norm_chrom_median, color = lineage)) +
        geom_hline(yintercept = c(1, 2), color = "black", linetype = "solid") +
        geom_point()+
        scale_x_continuous(limits = c(0,100))+
        theme_bw() +
        theme(legend.position = "none", legend.direction = "vertical")+
        labs(y = "Normalized Chromosome Median",
            x = "Percent of Chromosome Covered by Duplications",
            color = "Lineage",
            title = "Normalized Median Depth of Chromosomes vs Percentage of Chromosome Covered by CNVs")

del <- ggplot(dels, aes(x = coverage_percent, y = norm_chrom_median, color = lineage)) +
        geom_hline(yintercept = c(1, 2), color = "black", linetype = "solid") +
        geom_point()+
        scale_x_continuous(limits = c(0,100))+
        theme_bw() +
        theme(legend.position = "bottom", legend.direction = "horizontal")+
        labs(y = "Normalized Chromosome Median",
            x = "Percent of Chromosome Covered by Deletions",
            color = "Lineage")

p <- dup / del

print("Saving plot...")
ggsave(snakemake@output[[1]], plot = p, units = "in", height = 9, width = 16, scale = gscale)
print("Done!")
