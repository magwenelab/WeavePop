log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(RColorBrewer)
suppressPackageStartupMessages(library(scales))
library(svglite)
metadata <- read.csv(snakemake@input[[2]], header = TRUE)%>%
    select(sample, strain, lineage = group)

# metadata <- read.csv("./sample_metadata.csv", header = TRUE)%>%
#    select(sample, strain, lineage = Group)

stats <- read.delim(snakemake@input[[1]], sep =":", header = FALSE, col.names = c("stat", "value", "sample"))

#stats <- read.delim("./results/mapping_stats.txt", sep =":", header = FALSE, col.names = c("stat", "value", "sample"))
stats <- stats %>% pivot_wider(names_from = stat, values_from = value)
colnames(stats) <- gsub(" ", "_", colnames(stats))
stats <- stats %>%
    mutate(total_reads = reads_mapped + reads_unmapped,
            percent_mapped = (reads_mapped/total_reads)*100)%>%
            as.data.frame()

stats <- left_join(stats, metadata, by="sample")

stats_long <- stats %>%
    select(sample, lineage, strain, mapped = percent_mapped, properly_paired = "percentage_of_properly_paired_reads_(%)") %>%
    pivot_longer(c(mapped, properly_paired), names_to = "measurement", values_to = "value")%>%
    mutate(name = paste(strain, sample, sep=" " ))


plot <- ggplot(stats_long, aes(color = measurement, x=name, y= value))+
    geom_point()+
    ylim(0,100)+
    facet_grid(~lineage, scale = "free_x" , space='free_x')+
    scale_color_discrete(labels=c('Mapped', 'Properly paired'),guide = guide_legend(title = NULL))+
    labs(shape = NULL)+
    theme_light()+
    xlab("Sample") +
    ylab("Percentage of reads") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5))

#ggsave("./results/mapped_reads.svg", plot = plot, dpi = 50, units = "cm", height = 15, width = 60)
ggsave(snakemake@output[[1]], plot = plot, dpi = 50, units = "cm", height = 15, width = 60)