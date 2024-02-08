log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))

# Reading files

raw<- read.delim(snakemake@input[[1]], header = FALSE, col.names = c("Accession", "Start", "End", "Depth"), stringsAsFactors = TRUE)
good<- read.delim(snakemake@input[[2]], header = FALSE, col.names = c("Accession", "Start", "End", "Depth"), stringsAsFactors = TRUE)
chrom_names <- read.csv(snakemake@input[[3]], header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))

# Adding Chromosoma names to good and raw
raw <- left_join(raw, chrom_names, by = "Accession")
good <- left_join(good, chrom_names, by = "Accession")

# Calculating global coverage, normalized coverage by position (window), and coverage by chromosome
good_stats_regions <- good %>%
  mutate(Global_Mean = round(mean(Depth),2),
         Global_Median = round(median(Depth),2))%>%
  mutate(Norm_Mean= round(Depth/Global_Mean, 2))%>%
  mutate(Norm_Median= round(Depth/Global_Median, 2))%>%
  group_by(Chromosome)%>%
  mutate(Chrom_Mean = round(mean(Depth),2),
         Chrom_Median = round(median(Depth),2))%>%
  ungroup()%>%
  select(Chromosome, Accession, Lineage, Start, End, Depth, Global_Mean, Global_Median, Chrom_Mean, Chrom_Median, Norm_Mean, Norm_Median)

write.table(good_stats_regions, snakemake@output[[1]], col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)

raw_stats_regions <- raw %>%
  mutate(Global_Mean = round(mean(Depth),2),
         Global_Median = round(median(Depth),2))%>%
  mutate(Norm_Mean= round(Depth/Global_Mean, 2))%>%
  mutate(Norm_Median= round(Depth/Global_Median, 2))%>%
  group_by(Chromosome)%>%
  mutate(Chrom_Mean = round(mean(Depth),2),
         Chrom_Median = round(median(Depth),2))%>%
  ungroup()%>%
  select(Chromosome, Accession, Lineage, Start, End, Depth, Global_Mean, Global_Median, Chrom_Mean, Chrom_Median, Norm_Mean, Norm_Median)


write.table(raw_stats_regions, snakemake@output[[2]], col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)




