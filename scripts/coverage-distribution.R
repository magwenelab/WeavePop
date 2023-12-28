log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(RColorBrewer)
suppressPackageStartupMessages(library(scales))
library(svglite)

print("Reading files and joining data with chromosome names")
#setwd("./genomes-annotations/SRS404449/")
# cov<- read.csv("cov.csv", header = TRUE, stringsAsFactors = TRUE)
# chrom_names <- read.csv("../../chromosome_names.csv", header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))
# sample <- "/analysis/czirion/Crypto_Diversity_Pipeline/genomes-annotations/SRS404449/cov.csv"
sample <- snakemake@input[[1]]
Split <- str_split(sample, "/")
sample <- Split[[1]][length(Split[[1]])-1]

cov<- read.csv(snakemake@input[[1]], header = TRUE, stringsAsFactors = TRUE)
chrom_names <- read.csv(snakemake@input[[2]], header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))
cov <-cov %>%
  rename(Accession = Chromosome)
cov <- left_join(cov, chrom_names, by = "Accession")

print("Plotting Coverage distribution")
data_color = "aquamarine4"
lineage <- levels(as.factor(cov$Lineage))

plot <- ggplot(cov, aes(x=Coverage, y=Count))+
  geom_col(fill=data_color)+ 
  facet_wrap(~Chromosome,ncol = 2)+
  scale_y_log10(name = "Number of Sites", labels = comma)+
  scale_x_continuous(name = "Coverage (X) ", labels = comma, n.breaks = 10)+
  theme_light()+
  theme(legend.position="none")+
  labs(title = paste(lineage, sample,  sep = " "))

#ggsave("../../cov_distribution.svg", plot = plot, dpi = 50, units = "cm", height = 22, width = 22)

print("Saving plot")
ggsave(snakemake@output[[1]], plot = plot, dpi = 50, units = "cm", height = 22, width = 22)
print("Done!")
