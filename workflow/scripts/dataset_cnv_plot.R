log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggrepel))

print("Reading input parameters...")

metadata_path <- snakemake@input$metadata
chromosomes_path <- snakemake@input$chromosomes
cnv_path <- snakemake@input$cnv

output <- snakemake@output$plot
output2 <- snakemake@output$plot2

# metadata_path <- "/PhastData/czirion/WeavePop_workdir/Cneoformans_global/rhodes/results/02.Dataset/metadata.csv"
# chromosomes_path <- "/PhastData/czirion/WeavePop_workdir/Cneoformans_global/rhodes/results/02.Dataset/chromosomes.csv"

# cnv_path <- "/PhastData/czirion/WeavePop_workdir/Cneoformans_global/rhodes/results/02.Dataset/cnv/cnv_chromosomes.tsv"

# output2 <- "plot.png"

threshold_std = 0.3
threshold_smiley = 0.3

print("Reading files...")
metadata <- read.delim(metadata_path, sep = ",", header = TRUE, stringsAsFactors = TRUE)

if ("strain" %in% colnames(metadata)){
    metadata <- metadata %>%
        select(sample, strain, ref_genome) %>%
        mutate(name = paste(strain, sample, sep = " ")) 
    } else {
    metadata <- metadata %>%
        select(sample, ref_genome) %>%
        mutate(name = sample)}
    
metadata <- metadata %>%
    arrange(name)%>%
    mutate(name = factor(name, levels = unique(name)))

chromosomes <- read.delim(chromosomes_path, sep = ",", header = TRUE, stringsAsFactors = TRUE)
cnv <- read.delim(cnv_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
cnv <- left_join(cnv, metadata, by = "sample")%>%
    left_join(chromosomes, by = c("accession", "ref_genome"))%>%
    mutate(cnv = factor(str_replace(str_to_title(cnv), "_", " "),levels = c("Deletion","Duplication", "Single copy")))%>%
    group_by(sample)%>%
    mutate(std_norm_depth = median(chrom_std_norm_depth),
            median_slope = median(slope),
            genome_depth = first(genome_depth),
            flag_std = ifelse(std_norm_depth > threshold_std, "Non-uniform", "."),
            flag_smiley = ifelse(median_slope > threshold_smiley, "Smiley", "."),
            flag = case_when(flag_std == "Non-uniform" & flag_smiley == "Smiley" ~ "Non-uniform Smiley",
                            flag_std == "Non-uniform" & flag_smiley == "." ~ "Non-uniform",
                            flag_std == "." & flag_smiley == "Smiley" ~ "Smiley",
                            TRUE ~ NA))%>%
    ungroup()%>%
    arrange(name)%>%
    mutate(color = ifelse(is.na(flag), "black", "red"),
            colored_label = paste0("<span style='color:", color, "'>", name, "</span>"))%>%
    mutate(colored_label = factor(colored_label, levels = unique(colored_label)))

set2 <- brewer.pal(8, "Set2") 
#c <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
s_colors <- set2[c(6,5,8)]
names(s_colors) <- levels(cnv$cnv)
flag_colors <- c("Non-uniform" = "red", "Smiley" = "blue", "Non-uniform Smiley" = "purple", "Uniform" = "black")

plot <-ggplot(data = cnv) +
        geom_bar(aes(x = colored_label, y = coverage_percent, fill = cnv), stat = "identity") +
        scale_fill_manual(values = s_colors, name = NULL)+
        facet_grid(chromosome ~ ref_genome, scales = "free", space = "free_x") +
        theme(panel.background = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
            axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1, vjust = 0.5, size = 5))+
        labs(x = "", y = "% Coverage", title = "Percentage of Chromosome Covered by CNVs",
        subtitle = "Sample names with not uniform depth\nor smiley-face pattern are shown in red.")

cnv_sample <- cnv %>%
    select(sample, genome_depth, std_norm_depth, median_slope, flag)%>%
    distinct()%>%
    mutate(flag = ifelse(is.na(flag), "Uniform", flag))

a <- ggplot(cnv_sample)+
    geom_point(aes(x = genome_depth, y = std_norm_depth,
                    color = flag))+
    scale_color_manual(values = flag_colors, name = NULL, na.value = "black") +
    geom_text_repel(aes(x = genome_depth, y = std_norm_depth,
                    label = sample), size = 3)+
    theme_bw()+
    labs(title = "Uniformity of normalized depth vs. Genome-wide depth",
     x = "Genome depth",
     y = "Std dev of normalized depth")

b <- ggplot(cnv_sample)+
    geom_point(aes(x = genome_depth, y = median_slope,
                    color = flag))+
    scale_color_manual(values = flag_colors, name = NULL, na.value = "black") +   
    geom_text_repel(aes(x = genome_depth, y = median_slope,
                    label = sample), size = 3)+
    theme_bw()+
    theme(legend.position = "none")+
    labs(title = "Smiley face pattern vs. Genome-wide depth",
     x = "Genome depth",
     y = "Median of slope of normalized depth\nvs. relative distance to center of chromosome")

plot2 <- a / b 

print("Saving plot...")
n_samples <- nrow(metadata)
n_chroms <- length(unique(chromosomes$chromosome))

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

ggsave(output, plot = plot, units = "in", height = gheight, width = gwidth, scale = gscale, limitsize = FALSE)
ggsave(output2, plot = plot2, units = "in", height = 10, width = 10, limitsize = FALSE)

print("Done!")