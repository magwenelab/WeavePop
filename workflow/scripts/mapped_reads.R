log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scales))

metadata <- read.csv(snakemake@input[[2]], header = TRUE)%>%
    select(sample, strain, lineage = group)

# metadata <- read.csv("config/sample_metadata.csv", header = TRUE)%>%
#    select(sample, strain, lineage = group)

stats <- read.delim(snakemake@input[[1]], sep =":", header = FALSE, col.names = c("stat", "value", "sample"))

# stats <- read.delim("./results/dataset/mapping_stats.txt", sep =":", header = FALSE, col.names = c("stat", "value", "sample"))
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
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 17))+
    xlab("Sample") +
    ylab("Percentage of reads")

gwidth <- 13.3
gheight <- 7.5

ggsave(snakemake@output[[1]], plot = plot, scale = 0.5 , units = "in", height = gheight, width = gwidth)