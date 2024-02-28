log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scales))

#metadata <- read.csv("config/sample_metadata.csv", header = TRUE, stringsAsFactors = TRUE)

metadata <- read.csv(snakemake@input[[1]], header = TRUE, stringsAsFactors = TRUE)
metadata <- mutate(metadata, name = paste(strain, sample, sep = " "))

# chrom_names <- read.csv("config/chromosome_names.csv", header = FALSE, col.names = c("group", "Accession", "Chromosome"), stringsAsFactors = TRUE)
chrom_names <- read.csv(snakemake@input[[4]], header = FALSE, col.names = c("group", "Accession", "Chromosome"), stringsAsFactors = TRUE)
chrom_names <- chrom_names %>%
    select(Accession, Chromosome)

# good_stats <- read.delim("results/dataset/files/coverage_good.tsv", sep = "\t", header = TRUE, stringsAsFactors = TRUE)

good_stats <- read.delim(snakemake@input[[2]], sep = "\t", header = TRUE, stringsAsFactors = TRUE)
good_stats <- left_join(good_stats, metadata, by = "sample")
good_stats <- left_join(good_stats, chrom_names, by = "Accession")
good_stats <- good_stats %>%
    group_by(Accession, sample) %>%
    mutate(Norm_Mean = round(Chrom_Mean / Global_Mean, 2)) %>%
    mutate(Norm_Median = round(Chrom_Median / Global_Median, 2)) %>%
    ungroup() %>%
    as.data.frame()
# raw_stats <- read.delim("results/dataset/files/coverage_raw.tsv", sep = "\t", header = TRUE, stringsAsFactors = TRUE)

raw_stats <- read.delim(snakemake@input[[3]], sep = "\t", header = TRUE, stringsAsFactors = TRUE)
raw_stats <- left_join(raw_stats, metadata, by = "sample")
raw_stats <- left_join(raw_stats, chrom_names, by = "Accession")

raw_stats <- raw_stats %>%
    group_by(Accession, sample) %>%
    mutate(Norm_Mean = round(Chrom_Mean / Global_Mean, 2)) %>%
    mutate(Norm_Median = round(Chrom_Median / Global_Median, 2)) %>%
    ungroup() %>%
    as.data.frame()

gwidth <- 13.3
gheight <- 7.5
if (nlevels(good_stats$sample) <= 15) {
    gwidth <- 6
    gheight <- 7
    gscale <- 1
    gsize <- 5
} else if (nlevels(good_stats$sample) > 15 && nlevels(good_stats$sample) <= 60) {
    gscale <- 0.8
    gsize <- 8
} else if (nlevels(good_stats$sample) > 60 && nlevels(good_stats$sample) <= 150) {
    gscale <- 1
    gsize <- 6
} else if (nlevels(good_stats$sample) > 150 && nlevels(good_stats$sample) <= 400) {
    gscale <- 1.5
    gsize <- 3
} else {
    gscale <- 2
    gsize <- 2
}

# Median by Chromosome

toplim <- ceiling(max(good_stats$Norm_Median))
values <- seq(0, toplim, by = 1)
ylabel <- "Normalized coverage"

medianplot <- ggplot(good_stats, aes(x = reorder(name, -Global_Mean, sum), y = Norm_Median)) +
    geom_point(aes(color = get(snakemake@params[[1]]))) +
    ylim(0, toplim) +
    facet_grid(scale = "free_x", space = "free_x", rows = vars(Chromosome), cols = vars(group)) +
    scale_color_brewer(palette = "Set2", name = str_to_title(snakemake@params[[1]])) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = gsize),
        panel.grid.minor = element_blank()
    ) +
    #   strip.text = element_text(size = 15),
    #   legend.text = element_text(size = 15),
    #   legend.title = element_text(size = 17),
    #   plot.title = element_text(size = 20),
    #   axis.title = element_text(size = 17))+
    labs(
        title = "Normalized median coverage of chromosomes",
        x = "Sample",
        y = ylabel
    )


ggsave(snakemake@output[[1]], plot = medianplot, scale = gscale, units = "in", height = gheight, width = gwidth)

# Mean by Chromosome

toplim <- ceiling(max(good_stats$Norm_Mean))
values <- seq(0, toplim, by = 1)
ylabel <- "Normalized coverage"

meanplot <- ggplot(good_stats, aes(x = reorder(name, -Global_Mean, sum), y = Norm_Mean)) +
    geom_point(aes(color = get(snakemake@params[[1]]))) +
    ylim(0, toplim) +
    facet_grid(scale = "free_x", space = "free_x", rows = vars(Chromosome), cols = vars(group)) +
    scale_color_brewer(palette = "Set2", name = str_to_title(snakemake@params[[1]])) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = gsize),
        panel.grid.minor = element_blank()
    ) +
    #   strip.text = element_text(size = 15),
    #   legend.text = element_text(size = 15),
    #   legend.title = element_text(size = 17),
    #   plot.title = element_text(size = 20),
    #   axis.title = element_text(size = 17))+
    labs(
        title = "Normalized mean coverage of chromosomes",
        x = "Sample",
        y = ylabel
    )

ggsave(snakemake@output[[2]], plot = meanplot, scale = gscale, units = "in", height = gheight, width = gwidth)

# Global
topylim <- max(good_stats$Global_Mean) + max(good_stats$Global_Mean / 10)
raw_color = "#B3B3B3"
good_color = "#666666" 
color_quality = c("Good quality alignments" = good_color, "All alignments" = raw_color)
shape_stat <- c("Mean" = 16, "Median" = 15)

all <- bind_rows(good_stats %>% mutate(Quality = "Good quality alignments"), raw_stats %>% mutate(Quality = "All alignments"))
all$name <- reorder(all$name, -all$Global_Mean, sum)
all <- all %>%
    select(sample, name, Mean = Global_Mean, Median = Global_Median, Quality, group)%>%
    pivot_longer(
        cols = c(Mean, Median),
        names_to = "Measurement",
        values_to = "Value"
        )%>%
    as.data.frame()

all
topylim <- max(all$Value) + max(all$Value) / 10
g <- ggplot(all) +
    geom_point(aes(x = name, y = Value, shape = Measurement, color = Quality)) +
    ylim(0, topylim) +
    facet_grid(~group, scale = "free_x", space = "free_x") +
    scale_color_manual(name= "Alignment quality", values= color_quality)+
    scale_shape_manual(values = c(16,15), name = NULL, labels = c("Mean", "Median"))+
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = gsize),
        panel.grid.minor = element_blank()
    ) +
    # strip.text = element_text(size = 15),
    # legend.text = element_text(size = 15),
    # legend.title = element_text(size = 17),
    # plot.title = element_text(size = 20),
    # axis.title = element_text(size = 17))+
    labs(
        title = "Genome-wide coverage",
        x = "Sample",
        y = "Coverage (X)"
    )

ggsave(snakemake@output[[3]], plot = g, scale = gscale, units = "in", height = gheight, width = gwidth)
