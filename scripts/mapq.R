log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(RColorBrewer)
suppressPackageStartupMessages(library(scales))
library(svglite)
library(ggnewscale)

print("Reading BED file")

sample <- snakemake@input[[1]]
Split <- str_split(sample, "/")
sample <- Split[[1]][length(Split[[1]])-1]

raw<- read.delim(snakemake@input[[1]], header = FALSE, col.names = c("Accession", "Start", "End", "MAPQ"), stringsAsFactors = TRUE)
chrom_names <- read.csv(snakemake@input[[2]], header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))
raw <- left_join(raw, chrom_names, by = "Accession")

loci <- read.delim(snakemake@input[[3]], header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
loci <- loci %>% rename(Accession = seq_id)
loci <- left_join(loci, chrom_names, by = "Accession")

lineage <- levels(as.factor(raw$Lineage))
loci <- loci %>% filter(Lineage %in% lineage)

raw_color <- "bisque3"

print("Plotting chromosome depth")
plot <- ggplot()+
  geom_point(data = raw, aes(x= Start, y = MAPQ), size = 0.5 , color = raw_color)+
  geom_point(data = loci, aes(x= start, y = 1, color = Loci), size = 1, shape = 15)+
  facet_wrap(~Chromosome,ncol = 2, scales = "free_x")+
  scale_x_continuous(name = "Position (bp) ", labels = comma)+
  theme_light()+
  theme(legend.position="right")+
  labs(title = paste(lineage, sample,  sep = " "), y = "Mapping quality (phred score)")

ggsave(snakemake@output[[1]], plot = plot, dpi = 50, units = "cm", height = 22, width = 22)

