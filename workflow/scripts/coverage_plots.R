log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(RColorBrewer))

sample <- snakemake@input[[1]]
Split <- str_split(sample, "/")
sample <- Split[[1]][length(Split[[1]])-1]

good_stats_regions <- read.delim(snakemake@input[[1]], sep= ",", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
raw_stats_regions <- read.delim(snakemake@input[[2]], sep= ",", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))

good_stats_chroms <- good_stats_regions %>%
  select(c(-Start, -End, -Depth, -Norm_Mean, -Norm_Median))%>%
  distinct()
good_stats_chroms$Sample <- sample
good_stats_chroms <- good_stats_chroms %>%
  select(Sample, Lineage, Chromosome, Global_Mean, Global_Median, Chrom_Mean, Chrom_Median)

write_csv(good_stats_chroms, snakemake@output[[3]], col_names = TRUE)

raw_stats_chroms <- raw_stats_regions %>%
  select(c(-Start, -End, -Depth, -Norm_Mean, -Norm_Median))%>%
  distinct()
raw_stats_chroms$Sample <- sample
raw_stats_chroms <- raw_stats_chroms %>%
  select(Sample, Lineage, Chromosome, Global_Mean, Global_Median, Chrom_Mean, Chrom_Median)

write_csv(good_stats_chroms, snakemake@output[[4]], col_names = TRUE)

good_stats_long <- good_stats_chroms %>%
  pivot_longer(c(Chrom_Mean, Chrom_Median), names_to = "Measurement", values_to = "Value")
raw_stats_long <- raw_stats_chroms %>%
  pivot_longer(c(Chrom_Mean, Chrom_Median), names_to = "Measurement", values_to = "Value")
toplim <- max(raw_stats_long$Value) + max(raw_stats_long$Value)/10

lineage <- levels(as.factor(good_stats_chroms$Lineage))

loci <- read.delim(snakemake@input[[3]], header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
loci <- loci %>% 
    select(Accession = seq_id, start, Loci)

chrom_names <- good_stats_regions %>%
    select(Accession, Chromosome, Lineage)%>%
    distinct()

print("Ploting good quality coverage")            
raw_color = "lightskyblue1"
good_color = "lightskyblue3"
color_quality = c("Good quality alignments" = good_color, "All alignments" = raw_color)
topCov <- quantile(good_stats_regions$Norm_Median, 0.75) * 4
good_stats_regions$Norm_Median <- ifelse(good_stats_regions$Norm_Median >= topCov, topCov, good_stats_regions$Norm_Median)
plot <- ggplot()+
  geom_col(data = good_stats_regions, aes(x= Start, y = Norm_Median), color = good_color)+ 
  facet_wrap(~Chromosome,ncol = 1, scales = "free_x")+
  scale_x_continuous(name = "Position (bp) ", labels = comma)+
  theme_bw()+
  theme(legend.position="right",panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  labs(title = paste(lineage, sample,  sep = " "),
        y = "Normalized coverage")

if (nrow(loci) != 0){
  loci_chrom <- left_join(loci, chrom_names, by = c("Accession"))
  loci_chrom <- loci_chrom %>% filter(Lineage %in% lineage)
  loci_colors <- colorRampPalette(brewer.pal(8, "Dark2")[-8])(nlevels(loci_chrom$Loci))
  plot <- plot +
  geom_point(data = loci_chrom, aes(x= start, y = 0, color = Loci), size = 3, shape = 15)+
  scale_color_manual(name = "Loci", values = loci_colors)
}

struc_vars <- read.delim(snakemake@input[[4]], sep= ",", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
struc_vars <- struc_vars %>%
  filter(!is.na(End))

if (nrow(struc_vars) != 0){
  structure_colors <- colorRampPalette(brewer.pal(11, "RdBu")[c(-1, -2)])(nlevels(struc_vars$Structure))
  
  plot <- plot +
    new_scale_color()+
    geom_segment(data = struc_vars, aes(x= Start, y = topCov -1, xend = End, yend = topCov -1, color = Structure), size = 3)+
    scale_color_manual(name = "Structural variant", values = structure_colors)
}

pheight <- 10 + 2 * length(unique((good_stats_regions$Chromosome)))
pwidth <- 25
ggsave(snakemake@output[[1]], plot = plot, units = "cm", height = pheight, width = pwidth)

plot <- ggplot()+
  ylim(0,toplim) +
  geom_hline(aes(yintercept = unique(raw_stats_long$Global_Median),linetype = "Global median", color= "All alignments"))+
  geom_hline(aes(yintercept = unique(raw_stats_long$Global_Mean), linetype = "Global mean", color= "All alignments"))+
  geom_hline(aes(yintercept = unique(good_stats_long$Global_Median),linetype = "Global median", color = "Good quality alignments"))+
  geom_hline(aes(yintercept = unique(good_stats_long$Global_Mean), linetype = "Global mean", color = "Good quality alignments"))+
  geom_point(data = raw_stats_long, aes(x = factor(Chromosome, levels = as.character(sort(unique(Chromosome)))), y = Value, shape = Measurement, color= "All alignments"))+ 
  geom_point(data = good_stats_long, aes(x = factor(Chromosome, levels = as.character(sort(unique(Chromosome)))), y = Value, shape = Measurement, color = "Good quality alignments"))+ 
  labs(y = "Coverage", x = "Chromosome", title = paste(lineage, sample,  sep = " "))+
  theme_bw()+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_shape_manual(values = c(16,15), name = NULL, labels = c("Mean", "Median"))+
  scale_linetype_manual(values = c("solid","dotted"), name = NULL)+
  scale_color_manual(name= "Alignment quality", values= color_quality)

ggsave(snakemake@output[[2]], plot = plot, units = "cm", height = 15, width = 15)
