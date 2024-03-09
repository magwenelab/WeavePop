log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(RColorBrewer)
suppressPackageStartupMessages(library(scales))

print("Reading files and joining data with chromosome names")
# setwd("./results/samples/samtools/SRS8318899")
# cov<- read.csv("distrib_cov.csv", header = TRUE, stringsAsFactors = TRUE)
# chrom_names <- read.csv("../../../../config/chromosome_names.csv", header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))
# sample <- "/mnt/FastData/czirion/DiversityPipeline/results/samples/samtools/SRS8318899/distrib_cov.csv"
sample <- snakemake@input[[1]]
Split <- str_split(sample, "/")
sample <- Split[[1]][length(Split[[1]])-1]

cov<- read.table(snakemake@input[[1]], header = TRUE, stringsAsFactors = TRUE, sep = "\t")
cov[is.na(cov)] <- 0
chrom_names <- read.csv(snakemake@input[[2]], header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))
cov <-cov %>%
  rename(Accession = Chromosome)
cov <- left_join(cov, chrom_names, by = "Accession")

print("Plotting Coverage distribution")
raw_color = "#B3B3B3"
good_color = "#666666" 
color_quality = c("Good quality alignments" = good_color, "All alignments" = raw_color)
lineage <- levels(as.factor(cov$Lineage))

plot <- ggplot(cov, aes(x=Coverage))+
  geom_col(aes(y = Count_Raw, fill= "All alignments"))+ 
  geom_col(aes(y = Count_Good, fill= "Good quality alignments"))+ 
  facet_wrap(~Chromosome,ncol = 2)+
  scale_y_log10(name = "Number of Sites", labels = comma)+
  scale_x_continuous(name = "Coverage (X) ", labels = comma, n.breaks = 10)+
  scale_fill_manual(name= "Alignment quality", values= color_quality)+
  theme(legend.position="none")+
  labs(title = paste(lineage, sample,  sep = " "))+
  theme_bw()

#ggsave("../../cov_distribution.png", plot = plot, units = "cm", height = 22, width = 22)

pheight <- 0.5 + length(unique((cov$Chromosome)))/2
pwidth <- pheight * 1.78
ggsave(snakemake@output[[1]], plot = plot, units = "in", height = pheight, width = pwidth, dpi = 600)

cov_global <- cov %>%
  select(Coverage, Count_Good, Count_Raw)%>%
  group_by(Coverage)%>%
  summarize(Count_Good_Global = sum(Count_Good), Count_Raw_Global = sum(Count_Raw))%>%
  ungroup()

plot <- ggplot(cov_global, aes(x=Coverage))+
  geom_col(aes(y = Count_Raw_Global, fill= "All alignments"))+ 
  geom_col(aes(y = Count_Good_Global, fill= "Good quality alignments"))+ 
  scale_y_log10(name = "Number of Sites", labels = comma)+
  scale_x_continuous(name = "Coverage (X) ", labels = comma, n.breaks = 10)+
  scale_fill_manual(name= "Alignment quality", values= color_quality)+
  theme(legend.position="none")+
  labs(title = paste(lineage, sample,  sep = " "))+
  theme_bw()

ggsave(snakemake@output[[2]], plot = plot,  units = "in", height = 4.5, width = 8, dpi = 600)
print("Done!")
