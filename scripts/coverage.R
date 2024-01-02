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

raw<- read.delim(snakemake@input[[1]], header = FALSE, col.names = c("Accession", "Start", "End", "Depth"), stringsAsFactors = TRUE)
good<- read.delim(snakemake@input[[2]], header = FALSE, col.names = c("Accession", "Start", "End", "Depth"), stringsAsFactors = TRUE)
chrom_names <- read.csv(snakemake@input[[3]], header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))
raw <- left_join(raw, chrom_names, by = "Accession")
good <- left_join(good, chrom_names, by = "Accession")

loci <- read.delim(snakemake@input[[4]], header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
loci <- loci %>% rename(Accession = seq_id)
loci <- left_join(loci, chrom_names, by = "Accession")

lineage <- levels(as.factor(raw$Lineage))
loci <- loci %>% filter(Lineage %in% lineage)

raw_color = "lightskyblue1"
good_color = "lightskyblue3"
color_quality = c("Good quality alignments" = good_color, "All alignments" = raw_color)

print("Plotting chromosome depth")
good <- good %>%
  filter(Depth != 0)%>%
  mutate(Global_Mean = round(mean(Depth),2),
         Global_Median = round(median(Depth),2))%>%
  mutate(pmean= round(Depth/Global_Mean, 2))%>%
  mutate(pmedian= round(Depth/Global_Median, 2))%>%
  ungroup()

topCov <- quantile(good$pmedian, 0.75) * 4
good$pmedian <- ifelse(good$pmedian >= topCov, topCov, good$pmedian)
plot <- ggplot()+
  geom_col(data = good, aes(x= Start, y = pmedian), color = good_color)+ 
  geom_point(data = loci, aes(x= start, y = 0, color = Loci), size = 1, shape = 15)+
  facet_wrap(~Chromosome,ncol = 2, scales = "free_x")+
  scale_x_continuous(name = "Position (bp) ", labels = comma)+
  theme_bw()+
  theme(legend.position="right",panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  labs(title = paste(lineage, sample,  sep = " "),
        y = "Normalized coverage")

ggsave(snakemake@output[[1]], plot = plot, units = "cm", height = 22, width = 22)

print("Getting mean and median")
stats_raw <- raw %>%
  group_by(Chromosome)%>%
  summarise(Mean = round(mean(Depth),2),
            Median = round(median(Depth),2))%>%
  pivot_longer(c(Mean, Median), names_to = "Measurement", values_to = "Value")

global_stats_raw<- raw %>%
  summarise(Mean = round(mean(Depth),2),
            Median = round(median(Depth),2))

stats_good <- good %>%
  group_by(Chromosome)%>%
  summarise(Mean = round(mean(Depth),2),
            Median = round(median(Depth),2))%>%
  pivot_longer(c(Mean, Median), names_to = "Measurement", values_to = "Value")

global_stats_good<- good %>%
  summarise(Mean = round(mean(Depth),2),
            Median = round(median(Depth),2))

toplim <- max(stats_raw$Value) + max(stats_raw$Value)/10

plot <- ggplot()+
  ylim(0,toplim) +
  geom_hline(aes(yintercept = global_stats_raw$Median,linetype = "Global median", color= "All alignments"))+
  geom_hline(aes(yintercept = global_stats_raw$Mean, linetype = "Global mean", color= "All alignments"))+
  geom_hline(aes(yintercept = global_stats_good$Median,linetype = "Global median", color = "Good quality alignments"))+
  geom_hline(aes(yintercept = global_stats_good$Mean, linetype = "Global mean", color = "Good quality alignments"))+
  geom_point(data = stats_raw, aes(x = factor(Chromosome, levels = as.character(sort(unique(Chromosome)))), y = Value, shape = Measurement, color= "All alignments"))+ 
  geom_point(data = stats_good, aes(x = factor(Chromosome, levels = as.character(sort(unique(Chromosome)))), y = Value, shape = Measurement, color = "Good quality alignments"))+ 
  labs(y = "Coverage", x = "Chromosome", title = paste(lineage, sample,  sep = " "))+
  theme_bw()+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_shape_manual(values = c(16,15), name = NULL)+
  scale_linetype_manual(values = c("solid","dotted"), name = NULL)+
  scale_color_manual(name= "Alignment quality", values= color_quality)

ggsave(snakemake@output[[2]], plot = plot, units = "cm", height = 15, width = 15)

stats_raw$Sample <- sample
stats_good$Sample <- sample
global_stats_raw$Sample <- sample
global_stats_good$Sample <- sample

write_csv(stats_raw, snakemake@output[[3]], col_names = FALSE)
write_csv(stats_good, snakemake@output[[4]], col_names = FALSE)
write_csv(global_stats_raw, snakemake@output[[5]], col_names = FALSE)
write_csv(global_stats_good, snakemake@output[[6]], col_names = FALSE)
print("Done!")