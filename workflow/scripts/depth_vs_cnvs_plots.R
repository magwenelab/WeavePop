log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))

print("Reading files...")
depth <- read.delim(snakemake@input[[1]], sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
cnv_calls <- read.delim(snakemake@input[[2]], sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
chrom_names <- read.csv(snakemake@input[[3]], sep = ",", header = FALSE, col.names = c("lineage", "accession", "chromosome"), stringsAsFactors = TRUE, na = c("", "N/A"))
chrom_lengths <- read.delim(snakemake@input[[4]], sep = "\t", header = FALSE, col.names = c("accession", "length"), stringsAsFactors = TRUE, na = c("", "N/A"))
sample <- snakemake@wildcards$sample

# depth_by_chrom_good_path <-
#     "test/results/01.Samples/depth_quality/sample1/depth_by_chrom_good.tsv"
# cnv_calls_path <-
#     "test/results/01.Samples/cnv/sample1/cnv_calls.tsv"
# cnv_calls <- read.delim(cnv_calls_path, sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
# depth <- read.delim(depth_by_chrom_good_path, sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
# chromosomes_path <-
#     "test/results/02.Dataset/chromosomes.csv"
# chrom_names <- read.csv(chromosomes_path, sep = ",", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
# chrom_lengths_path <-
#     "test/results/04.Intermediate_files/03.References/chromosome_lengths.tsv"
# chrom_lengths <- read.delim(chrom_lengths_path, sep = "\t", header = FALSE, col.names = c("accession", "length"), stringsAsFactors = TRUE, na = c("", "N/A"))

# sample <- "sample1"

print("Processing chromosome information...")
chrom_names$chromosome <- factor(chrom_names$chromosome, levels = unique(chrom_names$chromosome))

chromosomes <- left_join(chrom_names, chrom_lengths, by = "accession")


depth <- depth %>%
    select(sample, accession, norm_chrom_median)%>%
    left_join(chromosomes, by = "accession")

print("Calculating chromosome metrics of CNVs...")

if (nrow(cnv_calls) > 0) {
    cnv_calls <- cnv_calls%>%
        group_by(cnv, accession)%>%
        summarise(total_cnv_size = sum(region_size)) %>%
        ungroup()%>%
        as.data.frame()
    del_chrom <- cnv_calls %>%
        filter(cnv == "deletion")
    dup_chrom <- cnv_calls %>%
        filter(cnv == "duplication")

    if (nrow(del_chrom) > 0) {
        print("Joining chromosome metrics of deletions with depth...")
        del_metrics <- left_join(depth, del_chrom, by = c("accession"))%>%
                        mutate(cnv = "deletion", 
                            coverage_percent = round((total_cnv_size / length) * 100, 2))
    } else {
        print("No deletions found.")
        del_metrics <- depth %>%
                        mutate(cnv = "deletion",
                                coverage_percent = 0)
    }

    if (nrow(dup_chrom) > 0) {
        print("Joining chromosome metrics of duplications with depth...")
        dup_metrics <- left_join(depth, dup_chrom, by = c("accession"))%>%
                        mutate(cnv = "duplication",
                                coverage_percent = round((total_cnv_size / length) * 100, 2))

    } else {
        print("No duplications found.")

        dup_metrics <- depth %>%
            mutate(cnv = "duplication",
                    coverage_percent = 0)

    }
} else {
    print("No CNVs found.")
    print("Creating empty CNV metrics...")
    del_metrics <- depth %>%
        mutate(cnv = "deletion",
                coverage_percent = 0)
    dup_metrics <- depth %>%
        mutate(cnv = "duplication",
                coverage_percent = 0)
}

print("Joining chromosome metrics of deletions and duplications...")
chrom_metrics <- bind_rows(del_metrics, dup_metrics)%>%
    mutate(coverage_percent = ifelse(is.na(coverage_percent), 0, coverage_percent))
        
lineage <- unique(chrom_metrics$lineage)

print("Plotting...")
p <- ggplot(chrom_metrics, aes(x = coverage_percent, y = norm_chrom_median, color = chromosome, shape = cnv)) +
        geom_hline(yintercept = c(0, 1, 2), color = "black", linetype = "solid") +
        geom_point() +
        geom_text_repel(aes(label = chromosome), size = 3, max.overlaps = 10, show.legend = FALSE) +
        scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
        theme_bw() +
        theme(legend.position = "right") +
        labs(title = "Normalized Depth vs. Percent of CNV Coverage per Chromosome",
            subtitle = paste("Lineage:", lineage, " Sample:", sample, sep = " "),
             y = "Normalized Median Depth of Chromosome",
             x = "Percent of Chromosome Covered by CNVs",
             color = "Chromosome",
             shape = "CNV")

print("Saving plot...")
ggsave(snakemake@output[[1]], p,  width = 8, height = 6, dpi = 300)
print("Done.")