log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggnewscale))
print("Reading files")


sample <- snakemake@input[[1]]
Split <- str_split(sample, "/")
sample <- Split[[1]][length(Split[[1]])-1]
# sample <- "SRS8318899"
# raw<- read.delim("results/samples/samtools/SRS8318899/mapq_window.bed", header = FALSE, col.names = c("Accession", "Start", "End", "MAPQ"), stringsAsFactors = TRUE)
raw<- read.delim(snakemake@input[[1]], header = FALSE, col.names = c("Accession", "Start", "End", "MAPQ"), stringsAsFactors = TRUE)
# chrom_names <- read.csv("config/chromosome_names.csv", header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))
chrom_names <- read.csv(snakemake@input[[4]], header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))
raw <- left_join(raw, chrom_names, by = "Accession")
# loci <- read.delim("results/dataset/files/loci_to_plot.tsv", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
loci <- read.delim(snakemake@input[[5]], header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
loci <- loci %>% rename(Accession = seq_id)
lineage <- levels(as.factor(raw$Lineage))

struc_vars <- read.delim(snakemake@input[[2]], sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
# struc_vars <- read.delim("results/samples/mosdepth/SRS8318899/good_structural_variants.tsv", sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
struc_vars <- left_join(struc_vars, chrom_names, by = "Accession")
structure <- struc_vars %>%
  select(Chromosome, Start, End, Structure)
s_lim <- max(raw$MAPQ) + 40
set2 <- rev(brewer.pal(8, "Set2")[1:6])
s_colors <- set2[1:nlevels(structure$Structure)]

repeats_table <- read.delim(snakemake@input[[3]], sep= "\t", header = FALSE, col.names = c("Accession", "Start", "End", "Repeat_type"), stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
# repeats_table <- read.delim("/FastData/czirion/DiversityPipeline/results/references/VNI/repeats/VNI_repeats.bed", sep= "\t", header = FALSE, col.names = c("Accession", "Start", "End", "Repeat_type"), stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
repeats<- left_join(repeats_table, chrom_names, by = "Accession")%>%
  select(Chromosome, Start, End, Repeat_type)

repeats$Repeat_type <- ifelse(repeats$Repeat_type == "Simple_repeat", "Simple repeat", "Others")
repeats$Repeat_type <- factor(repeats$Repeat_type, levels = c("Simple repeat", "Others"))
r_lim <- max(raw$MAPQ) + 60
r_colors <- colorRampPalette(brewer.pal(12, "Paired"))(nlevels(repeats$Repeat_type))

raw_color <- "black"
palette1 <- brewer.pal(n = 7, name = "Dark2")
palette2 <- brewer.pal(n = 7, name = "Set2")
combined_palette <- c(palette1, palette2)

print("Plotting chromosome MAPQ")
plot <- ggplot()+
  coord_cartesian(ylim=c(0,r_lim + 20), xlim= c(0,max(raw$End)))+
  geom_point(data = raw, aes(x= Start, y = MAPQ), size = 0.5 , color = raw_color)+
  geom_segment(data = repeats, aes(x = Start, xend = End, y = r_lim, yend = r_lim, color = Repeat_type), linewidth = 2)+
    scale_color_manual(name = "Type of repetitive sequence", values = r_colors)+
    guides(color = guide_legend(order=1))+
    new_scale_color()+
  geom_segment(data = structure, aes(x = Start, xend = End, y = s_lim, yend = s_lim, color = Structure), linewidth = 2)+
    scale_color_manual(name = "Structural variant", values = s_colors)+
    guides(color = guide_legend(order=2))+
    new_scale_color()+
  facet_wrap(~Chromosome,ncol = 2,strip.position = "right")+
  scale_x_continuous(name = "Position (bp) ", labels = comma)+
  labs(title = paste("Lineage: ", lineage,"Sample: ", sample,  sep = " "), y = "Mapping quality (phred score)")+
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 2))+
  scale_y_continuous(breaks = c(0, 20, 40, 60))

y_loci = max(raw$MAPQ) + 20
if (nrow(loci) != 0){
  loci_chrom <- left_join(loci, chrom_names, by = "Accession")
  loci_chrom <- loci_chrom %>% filter(Lineage %in% lineage)
  loci_colors <- combined_palette[1:nlevels(loci_chrom$Loci)]
  plot <- plot +
    geom_point(data = loci_chrom, aes(x= start, y = y_loci, color = Loci), size = 2, shape = 19)+
    scale_color_manual(name = "Loci", values = loci_colors)
}

ggsave(snakemake@output[[1]], plot = plot, units = "in", height = 9, width = 16)

