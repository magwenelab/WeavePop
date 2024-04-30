log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggnewscale))

# sample <- "SRS8318899"
# raw<- read.delim("results/samples/samtools/SRS8318899/mapq_window.bed", header = FALSE, col.names = c("Accession", "Start", "End", "MAPQ"), stringsAsFactors = TRUE)
# chrom_names <- read.csv("config/chromosome_names.csv", header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))
# loci <- read.delim("results/dataset/files/loci_to_plot.tsv", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
# struc_vars <- read.delim("results/samples/mosdepth/SRS8318899/good_structural_variants.tsv", sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
# repeats_table <- read.delim("/FastData/czirion/DiversityPipeline/results/references/VNI/repeats/VNI_repeats.bed", sep= "\t", header = FALSE, col.names = c("Accession", "Start", "End", "Repeat_type"), stringsAsFactors = TRUE, na = c("", "N/A", "NA"))

print("Reading files")
raw<- read.delim(snakemake@input[[1]], header = FALSE, col.names = c("Accession", "Start", "End", "MAPQ"), stringsAsFactors = TRUE)
struc_vars <- read.delim(snakemake@input[[2]], sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
repeats_table <- read.delim(snakemake@input[[3]], sep= "\t", header = FALSE, col.names = c("Accession", "Start", "End", "Repeat_type"), stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
chrom_names <- read.csv(snakemake@input[[4]], header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))
loci <- read.delim(snakemake@input[[5]], header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))

sample <- snakemake@wildcards$sample

chrom_names <- chrom_names %>%
  filter(Accession %in% unique(raw$Accession) )
chrom_names['Accession_Chromosome'] <- paste(chrom_names$Chromosome, chrom_names$Accession, sep = "xxx")
unique_levels <- unique(chrom_names$Accession_Chromosome)
new_order <- c(rbind(matrix(unique_levels, nrow = 2, byrow = TRUE)))
chrom_names$Accession_Chromosome <- factor(chrom_names$Accession_Chromosome, levels = new_order)

raw <- left_join(raw, chrom_names, by = "Accession")
loci <- loci %>% rename(Accession = seq_id)
lineage <- levels(as.factor(raw$Lineage))

struc_vars <- left_join(struc_vars, chrom_names, by = "Accession")
structure <- struc_vars %>%
  select(Accession_Chromosome, Chromosome, Start, End, Structure)
s_lim <- max(raw$MAPQ) + 40
set2 <- rev(brewer.pal(8, "Set2")[1:6])
s_colors <- set2[1:nlevels(structure$Structure)]

repeats<- left_join(repeats_table, chrom_names, by = "Accession")%>%
  select(Accession_Chromosome, Chromosome, Start, End, Repeat_type)

repeats$Repeat_type <- ifelse(repeats$Repeat_type == "Simple_repeat", "Simple repeat", "Others")
repeats$Repeat_type <- factor(repeats$Repeat_type, levels = c("Simple repeat", "Others"))
r_lim <- max(raw$MAPQ) + 60
r_colors <- colorRampPalette(brewer.pal(12, "Paired"))(nlevels(repeats$Repeat_type))

raw_color <- "black"
palette1 <- brewer.pal(n = 7, name = "Dark2")
palette2 <- brewer.pal(n = 7, name = "Set2")
combined_palette <- c(palette1, palette2)


my_labeller <- function(value){
  new_value <- sapply(strsplit(as.character(value), "xxx"), head, 1)
  return(new_value)
}

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
  facet_wrap(~Accession_Chromosome, strip.position = "right", ncol = 2, labeller = as_labeller(my_labeller)) +
  scale_x_continuous(name = "Position (bp) ", labels = comma)+
  labs(title = paste("Lineage: ", lineage," Sample: ", sample,  sep = " "), y = "Mapping quality (phred score)")+
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

ggsave(snakemake@output[[1]], plot = plot, units = "in", height = 9, width = 16, dpi = 600)

