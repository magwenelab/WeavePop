log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(RColorBrewer))


print("Reading input parameters...")

metadata_path <- snakemake@input$metadata
chromosomes_path <- snakemake@input$chromosomes
cnv_path <- snakemake@input$cnv

output <- snakemake@output$plot

# metadata_path <- "/FastData/czirion/WeavePop/test/results/02.Dataset/metadata.csv"
# chromosomes_path <- "/FastData/czirion/WeavePop/test/results/02.Dataset/chromosomes.csv"

# cnv_path <- "/FastData/czirion/WeavePop/test/results/02.Dataset/cnv/cnv_chromosomes.tsv"

# output <- "plot.png"

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

chromosomes <- read.delim(chromosomes_path, sep = ",", header = TRUE, stringsAsFactors = TRUE)
cnv <- read.delim(cnv_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
cnv <- left_join(cnv, metadata, by = "sample")%>%
    left_join(chromosomes, by = c("accession", "ref_genome"))%>%
    mutate(cnv = as.factor(str_replace(str_to_title(cnv), "_", " ")))
    
cnv$cnv <- factor(cnv$cnv, levels = c("Deletion","Duplication", "Single copy"))

set2 <- brewer.pal(8, "Set2") 
#c <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
s_colors <- set2[c(6,5,8)]
names(s_colors) <- levels(cnv$cnv)

plot <-ggplot() +
        geom_bar(data = cnv, aes(x = name, y = coverage_percent, fill = cnv), stat = "identity") +
        facet_grid(chromosome ~ ref_genome, scales = "free", space = "free_x") +
        theme(panel.background = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5))+
        labs(x = "", y = "% Coverage", title = "Percentage of Chromosome Covered by CNVs")+
        scale_fill_manual(values = s_colors, name = "")

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

print("Done!")