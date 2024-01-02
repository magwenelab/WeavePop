log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(RColorBrewer)
suppressPackageStartupMessages(library(scales))
library(svglite)

print("Reading files and joining data with chromosome names")
#setwd("./genomes-annotations/SRS404449/")
# mapq<- read.csv("mapq.csv", header = TRUE, stringsAsFactors = TRUE)
# chrom_names <- read.csv("../../chromosome_names.csv", header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))
# sample <- "/analysis/czirion/Crypto_Diversity_Pipeline/genomes-annotations/SRS404449/mapq.csv"
sample <- snakemake@input[[1]]
Split <- str_split(sample, "/")
sample <- Split[[1]][length(Split[[1]])-1]

mapq<- read.csv(snakemake@input[[1]], header = TRUE, stringsAsFactors = TRUE)
chrom_names <- read.csv(snakemake@input[[2]], header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))
mapq <-mapq %>%
  rename(Accession = Chromosome)
mapq <- left_join(mapq, chrom_names, by = "Accession")

print("Plotting MAPQ distribution")
data_color = "aquamarine4"
lineage <- levels(as.factor(mapq$Lineage))

plot <- ggplot(mapq, aes(x=MAPQ, y=Count))+
  geom_col(fill=data_color)+
  facet_wrap(~Chromosome,ncol = 2)+
  scale_y_log10(name = "Number of Reads", labels = comma, breaks = 100 * 10^seq(0,4, by = 2))+
  scale_x_continuous(name = "Mapping Quality", n.breaks = 8)+
  theme_bw()+
  theme(legend.position="none")+
  labs(title = paste(lineage, sample,  sep = " "))

#ggsave("../../mapq_distribution.svg", plot = plot, units = "cm", height = 22, width = 22)

print("Saving plot")
ggsave(snakemake@output[[1]], plot = plot, units = "cm", height = 22, width = 22)
print("Done!")

# Binned version

#mapq_bin <- mapq %>% 
#  mutate(bin_mapq = ntile(MAPQ,n=6))%>%
#  group_by(Chromosome, bin_mapq)%>%
#  reframe(binned_count = sum(Count))

#plot <- ggplot(mapq_bin, aes(x=bin_mapq, y=binned_count))+
#  geom_col(aes(fill = Chromosome))+
#  #scale_fill_manual(values = chrom_colors)+
#  facet_wrap(~Chromosome,ncol = 2)+
#  scale_y_log10(name = "Number of Reads", labels = comma, breaks = 100 * 10^seq(0,4, by = 2))+
#  scale_x_continuous(name = "Mapping Quality", n.breaks = 8)+
#  theme_bw()+
#  theme(legend.position="none")

##ggsave("mapq_distribution_binned.svg", plot = plot, units = "cm", height = 22, width = 22)