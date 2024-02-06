log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(RColorBrewer)
suppressPackageStartupMessages(library(scales))

print("Reading files")

sample <- snakemake@input[[1]]
Split <- str_split(sample, "/")
sample <- Split[[1]][length(Split[[1]])-1]

raw<- read.delim(snakemake@input[[1]], header = FALSE, col.names = c("Accession", "Start", "End", "MAPQ"), stringsAsFactors = TRUE)
chrom_names <- read.csv(snakemake@input[[2]], header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))
raw <- left_join(raw, chrom_names, by = "Accession")

loci <- read.delim(snakemake@input[[3]], header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
loci <- loci %>% rename(Accession = seq_id)
lineage <- levels(as.factor(raw$Lineage))

raw_color <- "#B3B3B3"
palette1 <- brewer.pal(n = 7, name = "Dark2")
palette2 <- brewer.pal(n = 7, name = "Set2")
combined_palette <- c(palette1, palette2)

print("Plotting chromosome MAPQ")
plot <- ggplot()+
  geom_point(data = raw, aes(x= Start, y = MAPQ), size = 0.5 , color = raw_color)+
  facet_wrap(~Chromosome,ncol = 2, scales = "free_x")+
  scale_x_continuous(name = "Position (bp) ", labels = comma)+
  theme_bw()+
  theme(legend.position="right")+
  labs(title = paste(lineage, sample,  sep = " "), y = "Mapping quality (phred score)")

if (nrow(loci) != 0){
  loci_chrom <- left_join(loci, chrom_names, by = "Accession")
  loci_chrom <- loci_chrom %>% filter(Lineage %in% lineage)
  loci_colors <- combined_palette[1:nlevels(loci_chrom$Loci)]
  plot <- plot +
    geom_point(data = loci_chrom, aes(x= start, y = 1, color = Loci), size = 2, shape = 19)+
    scale_color_manual(name = "Loci", values = loci_colors)
}

pheight <- 0.5 + length(unique((raw$Chromosome)))/2
pwidth <- pheight * 1.78
ggsave(snakemake@output[[1]], plot = plot, units = "in", height = pheight, width = pwidth)

