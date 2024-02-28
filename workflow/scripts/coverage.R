log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))

# Reading files
sample <- snakemake@input[[1]]
Split <- str_split(sample, "/")
sample <- Split[[1]][length(Split[[1]])-1]

raw<- read.delim(snakemake@input[[1]], header = FALSE, col.names = c("Accession", "Start", "End", "Depth"), stringsAsFactors = TRUE)
good<- read.delim(snakemake@input[[2]], header = FALSE, col.names = c("Accession", "Start", "End", "Depth"), stringsAsFactors = TRUE)

# Calculating global coverage, normalized coverage by position (window), and coverage by chromosome
good_stats_regions <- good %>%
  mutate(Global_Mean = round(mean(Depth),2),
         Global_Median = round(median(Depth),2))%>%
  mutate(Norm_Mean= round(Depth/Global_Mean, 2))%>%
  mutate(Norm_Median= round(Depth/Global_Median, 2))%>%
  group_by(Accession)%>%
  mutate(Chrom_Mean = round(mean(Depth),2),
         Chrom_Median = round(median(Depth),2))%>%
  ungroup()%>%
  select(Accession, Start, End, Depth, Global_Mean, Global_Median, Chrom_Mean, Chrom_Median, Norm_Mean, Norm_Median)

write.table(good_stats_regions, snakemake@output[[1]], col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

raw_stats_regions <- raw %>%
  mutate(Global_Mean = round(mean(Depth),2),
         Global_Median = round(median(Depth),2))%>%
  mutate(Norm_Mean= round(Depth/Global_Mean, 2))%>%
  mutate(Norm_Median= round(Depth/Global_Median, 2))%>%
  group_by(Accession)%>%
  mutate(Chrom_Mean = round(mean(Depth),2),
         Chrom_Median = round(median(Depth),2))%>%
  ungroup()%>%
  select(Accession, Start, End, Depth, Global_Mean, Global_Median, Chrom_Mean, Chrom_Median, Norm_Mean, Norm_Median)


write.table(raw_stats_regions, snakemake@output[[2]], col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

 # Getting stats per chromosome
good_stats_chroms <- good_stats_regions %>%
  select(c(-Start, -End, -Depth, -Norm_Mean, -Norm_Median))%>%
  distinct()
good_stats_chroms$sample <- sample
good_stats_chroms <- good_stats_chroms %>%
  select(sample, Accession, Global_Mean, Global_Median, Chrom_Mean, Chrom_Median)

write.table(good_stats_chroms, snakemake@output[[3]], col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

raw_stats_chroms <- raw_stats_regions %>%
  select(c(-Start, -End, -Depth, -Norm_Mean, -Norm_Median))%>%
  distinct()
raw_stats_chroms$sample <- sample
raw_stats_chroms <- raw_stats_chroms %>%
  select(sample, Accession, Global_Mean, Global_Median, Chrom_Mean, Chrom_Median)

write.table(raw_stats_chroms, snakemake@output[[4]], col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


