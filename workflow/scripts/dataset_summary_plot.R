log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(patchwork))

#metadata <- read.csv("config/sample_metadata.csv", header = TRUE, stringsAsFactors = TRUE)
metadata <- read.csv(snakemake@input[[1]], header = TRUE, stringsAsFactors = TRUE)
metadata <- mutate(metadata, name = paste(strain, sample, sep = " "))

# chrom_names <- read.csv("config/chromosome_names.csv", header = FALSE, col.names = c("lineage", "Accession", "Chromosome"), stringsAsFactors = TRUE)
chrom_names <- read.csv(snakemake@input[[2]], header = FALSE, col.names = c("lineage", "Accession", "Chromosome"), stringsAsFactors = TRUE)
chrom_names <- chrom_names %>%
    select(Accession, Chromosome)

# good_stats <- read.delim("results/dataset/files/coverage_good.tsv", sep = "\t", header = TRUE, stringsAsFactors = TRUE)
good_stats <- read.delim(snakemake@input[[3]], sep = "\t", header = TRUE, stringsAsFactors = TRUE)
good_stats <- rename(good_stats, sample = Sample)
good_stats <- left_join(good_stats, metadata, by = "sample")
good_stats <- left_join(good_stats, chrom_names, by = "Accession")
good_stats <- good_stats %>%
                group_by(Accession, sample) %>%
                mutate(Norm_Mean = round(Chrom_Mean / Global_Mean, 2)) %>%
                mutate(Norm_Median = round(Chrom_Median / Global_Median, 2)) %>%
                ungroup() %>%
                as.data.frame()

# raw_stats <- read.delim("results/dataset/files/coverage_raw.tsv", sep = "\t", header = TRUE, stringsAsFactors = TRUE)
raw_stats <- read.delim(snakemake@input[[4]], sep = "\t", header = TRUE, stringsAsFactors = TRUE)
raw_stats <- rename(raw_stats, sample = Sample)
raw_stats <- left_join(raw_stats, metadata, by = "sample")
raw_stats <- left_join(raw_stats, chrom_names, by = "Accession")

raw_stats <- raw_stats %>%
                group_by(Accession, sample) %>%
                mutate(Norm_Mean = round(Chrom_Mean / Global_Mean, 2)) %>%
                mutate(Norm_Median = round(Chrom_Median / Global_Median, 2)) %>%
                ungroup() %>%
                as.data.frame()
topylim <- max(good_stats$Global_Mean) + max(good_stats$Global_Mean / 10)
raw_color = "#B3B3B3"
good_color = "#666666" 
color_quality = c("Good quality mappings" = good_color, "All mappings" = raw_color)
shape_stat <- c("Mean" = 16, "Median" = 15)

all <- bind_rows(good_stats %>% mutate(Quality = "Good quality mappings"), raw_stats %>% mutate(Quality = "All mappings"))
all$name <- reorder(all$name, -all$Global_Mean, sum)
all <- all %>%
        select(sample, name, Mean = Global_Mean, Median = Global_Median, Quality, lineage)%>%
        pivot_longer(cols = c(Mean, Median), names_to = "Measurement", values_to = "Value")%>%
        as.data.frame()

topylim <- max(all$Value) + max(all$Value) / 10
g <- ggplot(all) +
    geom_point(aes(x = name, y = Value, shape = Measurement, color = Quality)) +
    ylim(0, topylim) +
    facet_grid(~lineage, scale = "free_x", space = "free_x") +
    scale_color_manual(name= "", values= color_quality)+
    scale_shape_manual(values = c(16,15), name = NULL, labels = c("Mean", "Median"))+
    theme_bw() +
    theme(panel.background = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    labs(title = "Genome-wide coverage",
         y = "Coverage (X)",
         x = "")

# Mapping stats and MAPQ
# map_stats <- read.table("results/dataset/files/mapping_stats.tsv",header = TRUE, stringsAsFactors = TRUE, sep = "\t")
map_stats <- read.table(snakemake@input[[5]], header = TRUE, stringsAsFactors = TRUE, sep = "\t")

metadata <- metadata %>%
    select(sample, strain, lineage)
stats_metad <- merge(map_stats, metadata, by = "sample")
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

name_order <- all %>%
    select(name)%>%
    distinct()%>%
    droplevels()
stats_reads$name <- factor(stats_reads$name, levels = levels(name_order$name))

stats_qualit <- stats_long %>%
    filter(Metric %in% c("percent_20", "percent_20_59", "percent_60"))
stats_qualit$Metric <- factor(stats_qualit$Metric, levels = c("percent_20","percent_20_59","percent_60"),
                              labels = c("MAPQ < 20", "MAPQ 20 - 59", "MAPQ 60"))

stats_qualit$name <- factor(stats_qualit$name, levels = levels(name_order$name))

palette_reads <- brewer.pal(n = length(unique(stats_reads$Metric)), name = "BuPu")
palette_qualit <- brewer.pal(n = length(unique(stats_qualit$Metric)), name = "BuGn")

reads <- ggplot()+
    geom_bar(data = stats_reads, aes(x = name, y = Value, fill = Metric), stat = "identity")+
    facet_grid(~ lineage, scales = "free", space = "free_x")+
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    labs(x = "", y = "Percentage of reads", fill = "Metric", title = "Percentage of reads by mapping status")+
    scale_fill_manual(values = palette_reads, name = "")

mapq <- ggplot()+
    geom_bar(data = stats_qualit, aes(x = name, y = Value, fill = Metric), stat = "identity")+
    facet_grid(~ lineage, scales = "free", space = "free_x")+
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = 5))+
    labs(x = "Sample", y = "Percentage of reads", fill = "Metric", title = "Percentage of mapped reads by mapping quality")+
    scale_fill_manual(values = palette_qualit, name = "")

plot <- g/reads/mapq

gscale = snakemake@params[[1]]
ggsave(snakemake@output[[1]], plot = plot, units = "in", height = 9, width = 16, scale = gscale)