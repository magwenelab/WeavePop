log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(RColorBrewer)
suppressPackageStartupMessages(library(scales))

metadata <- read.csv(snakemake@input[[1]], header = TRUE, stringsAsFactors = TRUE)
metadata <- mutate(metadata, name = paste(strain, sample, sep=" " ))

#### Good quality mappings ####
good_stats <-read.csv(snakemake@input[[2]], header = FALSE, col.names = c("sample", "Lineage", "Chromosome", "Global_Mean", "Global_Median", "Mean", "Median"), stringsAsFactors = TRUE)
good_stats <- left_join(good_stats, metadata, by = "sample")
good_stats <- good_stats%>%
    group_by(Chromosome, sample)%>%
    mutate(Norm_Mean= round(Mean/Global_Mean, 2))%>%
    mutate(Norm_Median= round(Median/Global_Median, 2))%>%
    ungroup()

write_csv(good_stats,snakemake@output[[1]], col_names = TRUE)

gwidth = 13.3
gheight = 7.5
if (nlevels(good_stats$sample) <= 15 ){
    gscale = 0.5
    gsize = 10
} else if (nlevels(good_stats$sample) > 15 && nlevels(good_stats$sample) <= 60 ){
    gscale = 0.8
    gsize = 8
} else if (nlevels(good_stats$sample) > 60 && nlevels(good_stats$sample) <= 150){
    gscale = 1 
    gsize = 6
} else if (nlevels(good_stats$sample) > 150 && nlevels(good_stats$sample)<= 400){
    gscale = 1.5
    gsize = 3
} else {
    gscale = 2
    gsize = 2
}

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
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = gsize),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 17))+
    labs(title= "Genome-wide coverage",
            x= "Sample",
            y= "Coverage (X)")



ggsave(snakemake@output[[2]], plot = g, scale = gscale,  units = "in", height = gheight, width = gwidth)

# Median by Chromosome 

toplim <- ceiling(max(good_stats$Norm_Median))
values <- seq(0, toplim, by = 1)
ylabel <- "Normalized coverage"

medianplot <- ggplot(good_stats, aes(x=reorder(name, -Global_Mean, sum), y= Norm_Median))+
    geom_point(aes(color= get(snakemake@params[[1]])))+
    ylim(0,toplim)+
    facet_grid(scale = "free_x" , space='free_x', rows= vars(Chromosome), cols = vars(group))+
    scale_color_brewer(palette = "Set2", name = str_to_title(snakemake@params[[1]]))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = gsize),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 17),
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 17))+
    labs(title = "Normalized median coverage of chromosomes",
         x = "Sample",
         y = ylabel)



ggsave(snakemake@output[[3]], plot = medianplot, scale = gscale,  units = "in", height = gheight, width = gwidth)

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
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = gsize),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 17),
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 17))+
    labs(title = "Normalized mean coverage of chromosomes",
         x = "Sample",
         y = ylabel)

ggsave(snakemake@output[[4]], plot = meanplot, scale = gscale,  units = "in", height = gheight, width = gwidth)

#### All quality mappings ####
raw_stats <-read.csv(snakemake@input[[2]], header = FALSE, col.names = c("sample", "Lineage", "Chromosome", "Global_Mean", "Global_Median", "Mean", "Median"), stringsAsFactors = TRUE)
raw_stats <- left_join(raw_stats, metadata, by = "sample")

raw_stats <- raw_stats%>%
    group_by(Chromosome, sample)%>%
    mutate(Norm_Mean= round(Mean/Global_Mean, 2))%>%
    mutate(Norm_Median= round(Median/Global_Median, 2))%>%
    ungroup()

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
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = gsize),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 17))+
    labs(title= "Genome-wide coverage",
            x= "Sample",
            y= "Coverage (X)")

ggsave(snakemake@output[[6]], plot = g, scale = gscale,  units = "in", height = gheight, width = gwidth)

# Median by Chromosome

toplim <- ceiling(max(raw_stats$Norm_Median))
values <- seq(0, toplim, by = 1)
ylabel <- "Normalized coverage"

medianplot <- ggplot(raw_stats, aes(x=reorder(name, -Global_Mean, sum), y= Norm_Median))+
    geom_point(aes(color= get(snakemake@params[[1]])))+
    ylim(0,toplim)+
    facet_grid(scale = "free_x" , space='free_x', rows= vars(Chromosome), cols = vars(group))+
    scale_color_brewer(palette = "Set2", name = str_to_title(snakemake@params[[1]]))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = gsize),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 17),
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 17))+
    labs(title = "Normalized median coverage of chromosomes",
         x = "Sample",
         y = ylabel)

ggsave(snakemake@output[[7]], plot = medianplot, scale = gscale,  units = "in", height = gheight, width = gwidth)

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
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = gsize),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 17),
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 17))+
    labs(title = "Normalized mean coverage of chromosomes",
         x = "Sample",
         y = ylabel)

ggsave(snakemake@output[[8]], plot = meanplot, scale = gscale,  units = "in", height = gheight, width = gwidth)