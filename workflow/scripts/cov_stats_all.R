log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(RColorBrewer)
suppressPackageStartupMessages(library(scales))

#metadata <- read.csv("files/sample_metadata.csv", header = TRUE, stringsAsFactors = TRUE)
metadata <- read.csv(snakemake@input[[1]], header = TRUE, stringsAsFactors = TRUE)
metadata <- mutate(metadata, name = paste(strain, sample, sep=" " ))

#### Good quality mappings ####
# good_stats <-read.csv("results/coverage_good.csv", header = FALSE, col.names = c("sample", "Lineage", "Chromosome", "Global_Mean", "Global_Median", "Mean", "Median"), stringsAsFactors = TRUE)
good_stats <-read.csv(snakemake@input[[2]], header = FALSE, col.names = c("sample", "Lineage", "Chromosome", "Global_Mean", "Global_Median", "Mean", "Median"), stringsAsFactors = TRUE)
good_stats <- left_join(good_stats, metadata, by = "sample")

good_stats <- good_stats%>%
    group_by(Chromosome, sample)%>%
    mutate(Norm_Mean= round(Mean/Global_Mean, 2))%>%
    mutate(Norm_Median= round(Median/Global_Median, 2))%>%
    ungroup()

#write_csv(chromosome,"results/norm_coverage_good.csv", col_names = TRUE)
write_csv(good_stats,snakemake@output[[1]], col_names = TRUE)


# Global
topylim <- max(good_stats$Global_Mean) + max(good_stats$Global_Mean/10)
color_stat = c("Mean" = "#008837", "Median" = "#7b3294")

g <- ggplot(good_stats, aes(x=reorder(name, -Global_Mean, sum)))+
    geom_point(aes(y= Global_Mean, color = "Mean"))+
    geom_point(aes(y= Global_Median, color = "Median"))+
    scale_color_manual(values= color_stat, name = "")+ 
    ylim(0,topylim)+
    facet_grid(~group,scale = "free_x" , space='free_x')+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 17))+
    labs(title= "Genome-wide coverage",
            x= "Sample",
            y= "Coverage (X)")

gwidth <- 12 + 0.15 * nlevels(good_stats$sample)
gheight <- gwidth/2

#ggsave("results/cov_global_good.png", plot = g,  units = "cm", height = gheight, width = gwidth)
ggsave(snakemake@output[[2]], plot = g,  units = "cm", height = gheight, width = gwidth)

# Median by Chromosome ####MAKE COLOR AES CONFIGURABLE ####

toplim <- ceiling(max(good_stats$Norm_Median))
values <- seq(0, toplim, by = 1)
ylabel <- "Normalized coverage"

medianplot <- ggplot(good_stats, aes(x=reorder(name, -Global_Mean, sum), y= Norm_Median))+
    geom_point(aes(color= get(snakemake@params[[1]])))+
    ylim(0,toplim)+
    facet_grid(scale = "free_x" , space='free_x', rows= vars(Chromosome), cols = vars(group))+
    scale_color_brewer(palette = "Set2", name = str_to_title(snakemake@params[[1]]))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 17),
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 17))+
    labs(title = "Normalized median coverage of chromosomes",
         x = "Sample",
         y = ylabel)

pwidth <- 12 + 0.15 * nlevels(good_stats$sample)
pheight <- 10 + 1 * length(unique(good_stats$Chromosome))

#ggsave("results/cov_median_good.png", plot = medianplot,  units = "cm", height = pheight, width = pwidth)
ggsave(snakemake@output[[3]], plot = medianplot,  units = "cm", height = pheight, width = pwidth)

# Mean by Chromosome

toplim <- ceiling(max(good_stats$Norm_Mean))
values <- seq(0, toplim, by = 1)
ylabel <- "Normalized coverage"

meanplot <- ggplot(good_stats, aes(x=reorder(name, -Global_Mean, sum), y= Norm_Mean))+
    geom_point(aes(color= get(snakemake@params[[1]])))+
    ylim(0,toplim)+
    facet_grid(scale = "free_x" , space='free_x', rows= vars(Chromosome), cols = vars(group))+
    scale_color_brewer(palette = "Set2", name = str_to_title(snakemake@params[[1]]))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 17),
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 17))+
    labs(title = "Normalized mean coverage of chromosomes",
         x = "Sample",
         y = ylabel)

#ggsave("results/cov_mean_good.png", plot = meanplot,  units = "cm", height = 30, width = 60)
ggsave(snakemake@output[[4]], plot = meanplot,  units = "cm", height = pheight, width = pwidth)

#### All quality mappings ####
# raw_stats <-read.csv("results/coverage_raw.csv", header = FALSE, col.names = c("sample", "Lineage", "Chromosome", "Global_Mean", "Global_Median", "Mean", "Median"), stringsAsFactors = TRUE)
raw_stats <-read.csv(snakemake@input[[2]], header = FALSE, col.names = c("sample", "Lineage", "Chromosome", "Global_Mean", "Global_Median", "Mean", "Median"), stringsAsFactors = TRUE)
raw_stats <- left_join(raw_stats, metadata, by = "sample")

raw_stats <- raw_stats%>%
    group_by(Chromosome, sample)%>%
    mutate(Norm_Mean= round(Mean/Global_Mean, 2))%>%
    mutate(Norm_Median= round(Median/Global_Median, 2))%>%
    ungroup()

#write_csv(chromosome,"results/norm_coverage_raw.csv", col_names = TRUE)
write_csv(raw_stats,snakemake@output[[5]], col_names = TRUE)


# Global
topylim <- max(raw_stats$Global_Mean) + max(raw_stats$Global_Mean/10)
color_stat = c("Mean" = "#008837", "Median" = "#7b3294")

g <- ggplot(raw_stats, aes(x=reorder(name, -Global_Mean, sum)))+
    geom_point(aes(y= Global_Mean, color = "Mean"))+
    geom_point(aes(y= Global_Median, color = "Median"))+
    scale_color_manual(values= color_stat, name = "")+ 
    ylim(0,topylim)+
    facet_grid(~group,scale = "free_x" , space='free_x')+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 17))+
    labs(title= "Genome-wide coverage",
            x= "Sample",
            y= "Coverage (X)")

#ggsave("results/cov_global_raw.png", plot = g,  units = "cm", height = gheight, width = gwidth)
ggsave(snakemake@output[[6]], plot = g,  units = "cm", height = gheight, width = gwidth)

# Median by Chromosome ####MAKE COLOR AES CONFIGURABLE ####

toplim <- ceiling(max(raw_stats$Norm_Median))
values <- seq(0, toplim, by = 1)
ylabel <- "Normalized coverage"

medianplot <- ggplot(raw_stats, aes(x=reorder(name, -Global_Mean, sum), y= Norm_Median))+
    geom_point(aes(color= get(snakemake@params[[1]])))+
    ylim(0,toplim)+
    facet_grid(scale = "free_x" , space='free_x', rows= vars(Chromosome), cols = vars(group))+
    scale_color_brewer(palette = "Set2", name = str_to_title(snakemake@params[[1]]))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 17),
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 17))+
    labs(title = "Normalized median coverage of chromosomes",
         x = "Sample",
         y = ylabel)

#ggsave("results/cov_median_raw.png", plot = medianplot,  units = "cm", height = pheight, width = pwidth)
ggsave(snakemake@output[[7]], plot = medianplot,  units = "cm", height = pheight, width = pwidth)

# Mean by Chromosome

toplim <- ceiling(max(raw_stats$Norm_Mean))
values <- seq(0, toplim, by = 1)
ylabel <- "Normalized coverage"

meanplot <- ggplot(raw_stats, aes(x=reorder(name, -Global_Mean, sum), y= Norm_Mean))+
    geom_point(aes(color= get(snakemake@params[[1]])))+
    ylim(0,toplim)+
    facet_grid(scale = "free_x" , space='free_x', rows= vars(Chromosome), cols = vars(group))+
    scale_color_brewer(palette = "Set2", name = str_to_title(snakemake@params[[1]]))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 17),
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 17))+
    labs(title = "Normalized mean coverage of chromosomes",
         x = "Sample",
         y = ylabel)

#ggsave("results/cov_mean_raw.png", plot = meanplot,  units = "cm", height = pheight, width = pwidth)
ggsave(snakemake@output[[8]], plot = meanplot,  units = "cm", height = pheight, width = pwidth)