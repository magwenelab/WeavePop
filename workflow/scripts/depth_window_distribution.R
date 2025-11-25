log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggbeeswarm))

print("Reading input parameters...")
windows_input <- snakemake@input$depth
chromosomes_input <- snakemake@input$chroms
metadata_input <- snakemake@input$metadata

output <- snakemake@output$plot

window_size <- snakemake@params$window_size
depth_threshold <- snakemake@params$depth_threshold
smoothing_size <- snakemake@params$smoothing_size
sample <- snakemake@wildcards$sample
gheight <- 9
gwidth <- 16
gdpi <- 600


print("Reading files...")
windows <- read.delim(windows_input, sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
chromosomes <- read.csv(chromosomes_input, sep= ",", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
metadata <- read.delim(metadata_input, sep = ",", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))

print("Obtaining lineage of sample...")

lineage_name <- as.character(metadata$lineage[metadata$sample == sample])
strain_name <- as.character(metadata$strain[metadata$sample == sample])

print("Obtaining lineage of sample...")

lineage_name <- as.character(metadata$lineage[metadata$sample == sample])

print("Filtering and ordering chromosomes...")
chromosomes <- chromosomes %>%
                filter(lineage == lineage_name)%>%
                arrange(length)

chromosomes$chromosome <- factor(chromosomes$chromosome, levels = unique(chromosomes$chromosome))

print("Joining data...")
windows <- left_join(windows, chromosomes, by = "accession")
windows$cnv <- factor(windows$cnv, levels = c("deletion", "duplication",  "single_copy"))

print("Plotting depth..")

set2 <- rev(brewer.pal(8, "Set2")[1:6])
s_colors <- c(set2[1:2], "black")

median_depth <- windows %>%
  summarise(median = median(depth)) %>%
  pull(median)

c <- ggplot(windows, aes(x = chromosome, y = depth))+
    geom_quasirandom(aes(color = cnv), alpha = 0.5)+
    scale_color_manual(values = s_colors, name = "CNV")+
    geom_boxplot(color = "gray", outlier.shape = NA, fill = NA)+
    scale_y_continuous(limits = c(0,median_depth*4))+
    theme_minimal()+
    theme(legend.position = "none",
        axis.title.x = element_blank())+
    labs(y="Depth\n(truncated)",
        title = "Distribution of all Depth Metrics of Windows per Chromosome and Whole Genome",
        subtitle = paste("Sample:", sample,
                    "Strain:", strain_name,
                    " Lineage:", lineage_name, 
                    " Window size:", window_size,  
                    " Depth threshold for CNV calling:", depth_threshold,
                    " Smoothing size for CNV calling:", smoothing_size,
                     sep = " "))

g <- ggplot(windows, aes(y = depth, x = 1))+
    geom_quasirandom(aes(color = cnv), alpha = 0.5)+
    scale_color_manual(values = s_colors, name = "CNV")+
    geom_boxplot(color = "gray", outlier.shape = NA, fill = NA)+
    scale_y_continuous(position = "right", limits = c(0,median_depth*4))+
    theme_minimal()+
    theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

d <- c + g + plot_layout(, nrow =1, widths = c(5,1))

print("Plotting normalized depth..")
c <- ggplot(windows, aes(x = chromosome, y = norm_depth))+
    geom_quasirandom(aes(color = cnv), alpha = 0.5)+
    scale_color_manual(values = s_colors, name = "CNV")+
    geom_boxplot(color = "gray", outlier.shape = NA, fill = NA)+
    scale_y_continuous( breaks = seq(0,5, by = 0.5), limits = c(0,4))+
    geom_hline(aes(yintercept = 1 + depth_threshold), color = "red", linetype = "dashed")+
    geom_hline(aes(yintercept = 1 - depth_threshold), color = "red", linetype = "dashed")+
    theme_minimal()+
    theme(legend.position = "none",
        axis.title.x = element_blank())+
    labs(y="Normalized Depth\n(truncated)")

g <- ggplot(windows, aes(y = norm_depth, x = 1))+
    geom_quasirandom(aes(color = cnv), alpha = 0.5)+
    scale_color_manual(values = s_colors, name = "CNV")+
    geom_boxplot(color = "gray", outlier.shape = NA, fill = NA)+
    scale_y_continuous( breaks = seq(0,5, by = 0.5), position = "right", limits = c(0,4))+
    geom_hline(aes(yintercept = 1 + depth_threshold), color = "red", linetype = "dashed")+
    geom_hline(aes(yintercept = 1 - depth_threshold), color = "red", linetype = "dashed")+
    theme_minimal()+
    theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")


n <- c + g + plot_layout(, nrow =1, widths = c(5,1))

print("Plotting smooth depth..")
c <- ggplot(windows, aes(x = chromosome, y = smooth_depth))+
    geom_quasirandom(aes(color = cnv), alpha = 0.5)+
    scale_color_manual(values = s_colors, name = "CNV")+
    geom_boxplot(color = "gray", outlier.shape = NA, fill = NA)+
    scale_y_continuous( breaks = seq(0,5, by = 0.5), limits = c(0,4))+
    geom_hline(aes(yintercept = 1 + depth_threshold), color = "red", linetype = "dashed")+
    geom_hline(aes(yintercept = 1 - depth_threshold), color = "red", linetype = "dashed")+
    theme_minimal()+
    theme(legend.position = "bottom")+
    labs(y="Smooth Normalized Depth\n(truncated)",
        x = "Chormosomes Ordered by Length")

g <- ggplot(windows, aes(y = smooth_depth, x = 1))+
    geom_quasirandom(aes(color = cnv), alpha = 0.5)+
    scale_color_manual(values = s_colors, name = "CNV")+
    geom_boxplot(color = "gray", outlier.shape = NA, fill = NA)+
    scale_y_continuous( breaks = seq(0,5, by = 0.5), position = "right", limits = c(0,4))+
    geom_hline(aes(yintercept = 1 + depth_threshold), color = "red", linetype = "dashed")+
    geom_hline(aes(yintercept = 1 - depth_threshold), color = "red", linetype = "dashed")+
    theme_minimal()+
    theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")+
    labs(x = "Whole Genome")

s <- c + g + plot_layout(, nrow =1, widths = c(5,1))

p <- d / n / s

print("Saving plot...")
ggsave(output, p, height = gheight, width = gwidth, dpi = gdpi, units = "in")

print("Done!")





