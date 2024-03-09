log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(patchwork))

metadata <- read.csv(snakemake@input[[2]], header = TRUE)%>%
    select(sample, strain, lineage = group)

# metadata <- read.csv("config/sample_metadata.csv", header = TRUE)%>%
#    select(sample, strain, lineage = group)
# stats <- read.table("results/dataset/files/mapping_stats.tsv",header = TRUE, stringsAsFactors = TRUE, sep = "\t")
stats <- read.table(snakemake@input[[1]],header = TRUE, stringsAsFactors = TRUE, sep = "\t")
stats_metad <- merge(stats, metadata, by = "sample")
stats_metad <- stats_metad %>%
    mutate(name = paste(strain,sample, sep = " "),
        reads_unmapped = raw_total_sequences - reads_mapped,
        percent_unmapped = (reads_unmapped/raw_total_sequences)*100,
        reads_only_mapped = reads_mapped - reads_properly_paired,
        percent_only_mapped = (reads_only_mapped/raw_total_sequences)*100,
        percent_properly_paired = (reads_properly_paired/raw_total_sequences)*100,
        percent_20 = (MAPQ_20/reads_mapped)*100,
        percent_20_59 = (MAPQ_20_59/reads_mapped)*100,
        percent_60 = (MAPQ_60/reads_mapped)*100)

stats_long <- stats_metad %>%
    pivot_longer(cols = -c(sample, name, lineage, strain), names_to = "Metric", values_to = "Value")
stats_reads <- stats_long %>%
    filter(Metric %in% c("percent_only_mapped", "percent_unmapped", "percent_properly_paired"))
stats_reads$Metric <- factor(stats_reads$Metric, levels = c("percent_unmapped", "percent_only_mapped", "percent_properly_paired"),
                              labels = c("Unmapped", "Mapped", "Mapped and properly paired"))
stats_qualit <- stats_long %>%
    filter(Metric %in% c("percent_20", "percent_20_59", "percent_60"))
stats_qualit$Metric <- factor(stats_qualit$Metric, levels = c("percent_20","percent_20_59","percent_60"),
                              labels = c("MAPQ < 20", "MAPQ 20 - 59", "MAPQ 60"))

palette_reads <- brewer.pal(n = length(unique(stats_reads$Metric)), name = "BuPu")
palette_qualit <- brewer.pal(n = length(unique(stats_qualit$Metric)), name = "BuGn")

reads <- ggplot()+
    geom_bar(data = stats_reads, aes(x = name, y = Value, fill = Metric), stat = "identity")+
    facet_grid(~ lineage, scales = "free", space = "free_x")+
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 2),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    labs(x = "", y = "Percentage of reads", fill = "Metric", title = "Percentage of reads by mapping classification")+
    scale_fill_manual(values = palette_reads, name = "")

mapq <- ggplot()+
    geom_bar(data = stats_qualit, aes(x = name, y = Value, fill = Metric), stat = "identity")+
    facet_grid(~ lineage, scales = "free", space = "free_x")+
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 2),
          axis.text.x = element_text(angle = 90))+
    labs(x = "", y = "Percentage of reads", fill = "Metric", title = "Percentage of mapped reads by low, medium and high mapping quality")+
    scale_fill_manual(values = palette_qualit, name = "")

plot <- reads / mapq 
ggsave(snakemake@output[[1]], plot = plot)
