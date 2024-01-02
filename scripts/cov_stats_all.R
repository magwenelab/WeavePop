log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(RColorBrewer)
suppressPackageStartupMessages(library(scales))
library(svglite)
library(ggnewscale)

#metadata <- read.csv("files/sample_metadata.csv", header = TRUE, stringsAsFactors = TRUE)
metadata <- read.csv(snakemake@input[[1]], header = TRUE, stringsAsFactors = TRUE)
metadata <- mutate(metadata, name = paste(strain, sample, sep=" " ))

#### Good quality mappings ####
# global <-read.csv("results/coverage_global_good.csv", header = FALSE, col.names = c("Global_Mean", "Global_Median", "sample"), stringsAsFactors = TRUE)
# chromosome <-read.csv("results/coverage_good.csv", header = FALSE, col.names = c("Chromosome", "Measurement", "Value", "sample"), stringsAsFactors = TRUE)
global <-read.csv(snakemake@input[[2]], header = FALSE, col.names = c("Global_Mean", "Global_Median", "sample"), stringsAsFactors = TRUE)
chromosome <-read.csv(snakemake@input[[3]], header = FALSE, col.names = c("Chromosome", "Measurement", "Value", "sample"), stringsAsFactors = TRUE)
global <- left_join(global, metadata, by = "sample")
chromosome$Chromosome <- as.factor(chromosome$Chromosome)

chromosome <- pivot_wider(chromosome, names_from = Measurement, values_from = Value)
chromosome <- left_join(chromosome, global, by = "sample") 
chromosome <- chromosome%>%
    group_by(Chromosome, sample)%>%
    mutate(pmean= round(Mean/Global_Mean, 2))%>%
    mutate(pmedian= round(Median/Global_Median, 2))%>%
    ungroup()

#write_csv(chromosome,"results/proportional_coverage_good.csv", col_names = TRUE)
write_csv(chromosome,snakemake@output[[1]], col_names = TRUE)


# Global
topylim <- max(global$Global_Mean) + max(global$Global_Mean/10)
color_stat = c("Mean" = "#008837", "Median" = "#7b3294")

g <- ggplot(global, aes(x=reorder(name, -Global_Mean, sum)))+
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

#ggsave("results/cov_good_all.svg", plot = g,  units = "cm", height = 30, width = 60)
ggsave(snakemake@output[[2]], plot = g,  units = "cm", height = 30, width = 60)

# Median by Chromosome

toplim <- ceiling(max(chromosome$pmedian))
values <- seq(0, toplim, by = 1)
#colors <- brewer.pal(n = toplim +1, name = "Dark2") 
ylabel <- "Normalized coverage"

medianplot <- ggplot(chromosome, aes(x=reorder(name, -Global_Mean, sum), y= pmedian))+
    geom_point(aes(color= source))+
    ylim(0,toplim)+
    facet_grid(scale = "free_x" , space='free_x', rows= vars(Chromosome), cols = vars(group))+
    #scale_color_gradientn(colors = colors, breaks = values,limits = c(0, toplim), values = rescale(values), guide = "colorbar", name = ylabel) +
    scale_color_brewer(palette = "Set2", name = "Source")+
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

#ggsave("results/cov_prop_median_good.svg", plot = medianplot,  units = "cm", height = 30, width = 60)
ggsave(snakemake@output[[3]], plot = medianplot,  units = "cm", height = 30, width = 60)

# Mean by Chromosome

toplim <- ceiling(max(chromosome$pmean))
values <- seq(0, toplim, by = 1)
#colors <- brewer.pal(n = toplim +1, name = "Dark2") 
ylabel <- "Normalized coverage"

meanplot <- ggplot(chromosome, aes(x=reorder(name, -Global_Mean, sum), y= pmean))+
    geom_point(aes(color= source))+
    ylim(0,toplim)+
    facet_grid(scale = "free_x" , space='free_x', rows= vars(Chromosome), cols = vars(group))+
    #scale_color_gradientn(colors = colors, breaks = values,limits = c(0, toplim), values = rescale(values), guide = "colorbar", name = ylabel) +
    scale_color_brewer(palette = "Set2", name = "Source")+
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

#ggsave("results/cov_prop_mean_good.svg", plot = meanplot,  units = "cm", height = 30, width = 60)
ggsave(snakemake@output[[4]], plot = meanplot,  units = "cm", height = 30, width = 60)

# lineage <- levels(chromosome$group)
# for (lin in lineage){ 
#     temp<- chromosome %>% filter(group == lin) %>% droplevels()
#     p <- ggplot(temp, aes(x=reorder(name, -Global_Mean, sum), y= pmedian))+
#         geom_point(aes(color= source), size =2)+
#         ylim(0,toplim)+
#         facet_grid(scale = "free_x" , space='free_x', rows= vars(Chromosome))+
#         scale_color_brewer(palette = "Set2", name = "Source")+
#         theme_bw()+
#         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6),
#             panel.grid.minor = element_blank(),
#             strip.text = element_text(size = 15),
#             legend.text = element_text(size = 15),
#             legend.title = element_text(size = 17),
#             plot.title = element_text(size = 20),
#             axis.title = element_text(size = 17))+
#         labs(title = paste("Normalized median coverage of chromosomes of lineage", lin, sep = " "),
#             x = "Sample",
#             y = ylabel)
#     assign(paste("p", lin, sep = ""), p)
# }

#     ggsave("results/cov_prop_median_good_VNI.svg", plot = pVNI,  units = "cm", height = 30, width = 45)
#     ggsave("results/cov_prop_median_good_VNII.svg", plot = pVNII,  units = "cm", height = 30, width = 15)
#     ggsave("results/cov_prop_median_good_VNBI.svg", plot = pVNBI,  units = "cm", height = 30, width = 30)
#     ggsave("results/cov_prop_median_good_VNBII.svg", plot = pVNBII,  units = "cm", height = 30, width = 27)

#### All mappings ####
# global <-read.csv("results/coverage_global_raw.csv", header = FALSE, col.names = c("Global_Mean", "Global_Median", "sample"), stringsAsFactors = TRUE)
# chromosome <-read.csv("results/coverage_raw.csv", header = FALSE, col.names = c("Chromosome", "Measurement", "Value", "sample"), stringsAsFactors = TRUE)
global <-read.csv(snakemake@input[[4]], header = FALSE, col.names = c("Global_Mean", "Global_Median", "sample"), stringsAsFactors = TRUE)
chromosome <-read.csv(snakemake@input[[5]], header = FALSE, col.names = c("Chromosome", "Measurement", "Value", "sample"), stringsAsFactors = TRUE)
global <- left_join(global, metadata, by = "sample")
chromosome$Chromosome <- as.factor(chromosome$Chromosome)

chromosome <- pivot_wider(chromosome, names_from = Measurement, values_from = Value)
chromosome <- left_join(chromosome, global, by = "sample") 
chromosome <- chromosome%>%
    group_by(Chromosome, sample)%>%
    mutate(pmean= round(Mean/Global_Mean, 2))%>%
    mutate(pmedian= round(Median/Global_Median, 2))%>%
    ungroup()

write_csv(chromosome,snakemake@output[[5]], col_names = TRUE)

# Global
topylim <- max(global$Global_Mean) + max(global$Global_Mean/10)
color_stat = c("Mean" = "#008837", "Median" = "#7b3294")

g <- ggplot(global, aes(x=reorder(name, -Global_Mean, sum)))+
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
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5))+
    labs(title= "Global coverage",
            x= "Sample",
            y= "Coverage (X)")

#ggsave("results/cov_raw_all.svg", plot = g,  units = "cm", height = 30, width = 60)
ggsave(snakemake@output[[6]], plot = g,  units = "cm", height = 30, width = 60)

# Median by Chromosome

toplim <- ceiling(max(chromosome$pmedian))
values <- seq(0, toplim, by = 1)
#colors <- brewer.pal(n = toplim +1, name = "Dark2") 
ylabel <- "Normalized coverage"

medianplot <- ggplot(chromosome, aes(x=reorder(name, -Global_Mean, sum), y= pmedian))+
    geom_point(aes(color= source))+
    ylim(0,toplim)+
    facet_grid(scale = "free_x" , space='free_x', rows= vars(Chromosome), cols = vars(group))+
    scale_color_brewer(palette = "Set2", name = "Source")+
    #scale_color_gradientn(colors = colors, breaks = values,limits = c(0, toplim), values = rescale(values), guide = "colorbar", name = ylabel) +
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

#ggsave("results/cov_prop_median_raw.svg", plot = medianplot,  units = "cm", height = 30, width = 60)
ggsave(snakemake@output[[7]], plot = medianplot,  units = "cm", height = 30, width = 60)

# Mean by Chromosome

toplim <- ceiling(max(chromosome$pmean))
values <- seq(0, toplim, by = 1)
#colors <- brewer.pal(n = toplim +1, name = "Dark2") 
ylabel <- "Normalized coverage"

meanplot <- ggplot(chromosome, aes(x=reorder(name, -Global_Mean, sum), y= pmean))+
    geom_point(aes(color= source))+
    ylim(0,toplim)+
    facet_grid(scale = "free_x" , space='free_x', rows= vars(Chromosome), cols = vars(group))+
    scale_color_brewer(palette = "Set2", name = "Source")+
    #scale_color_gradientn(colors = colors, breaks = values,limits = c(0, toplim), values = rescale(values), guide = "colorbar", name = ylabel) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 17),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 17))+
    labs(title = "Normalized mean coverage of chromosome",
         x = "Sample",
         y = ylabel)
#ggsave("results/cov_prop_mean_raw.svg", plot = medianplot,  units = "cm", height = 30, width = 60)

ggsave(snakemake@output[[8]], plot = meanplot,  units = "cm", height = 30, width = 60)