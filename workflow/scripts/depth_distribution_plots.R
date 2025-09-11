log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(patchwork))

print("Reading input parameters...")
depth_input <- snakemake@input$distrib
chrom_names_input <- snakemake@input$chroms
metadata_input <- snakemake@input$metadata

output_chrom <- snakemake@output$chrom
output_global <- snakemake@output$globl

sample <- snakemake@wildcards$sample
gscale <- 1.5
gheight <- 4.5
gwidth <- 8
gdpi <- 300


print("Reading files and joining data with chromosome names...")
depth <- read.table(depth_input, header = TRUE, stringsAsFactors = TRUE, sep = "\t")
chrom_names <- read.csv(chrom_names_input, header = TRUE , sep = ",")
metadata <- read.delim(metadata_input, sep = ",", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))

print("Obtaining lineage of sample...")

lineage_name <- as.character(metadata$lineage[metadata$sample == sample])
strain_name <- as.character(metadata$strain[metadata$sample == sample])

print("Adding 0 to null depth values...")
depth[is.na(depth)] <- 0

print("Filtering chromosome names...")
chrom_names <- chrom_names %>%
  filter(accession %in% unique(depth$accession) )

chrom_names$chromosome <- factor(chrom_names$chromosome, levels = unique(chrom_names$chromosome))

print("Joining and arranging data...")
depth <- left_join(depth, chrom_names, by = "accession")

depth_global <- depth %>%
  select(depth, count_good, count_raw)%>%
  group_by(depth)%>%
  summarize(count_good_global = sum(count_good), count_raw_global = sum(count_raw))%>%
  ungroup()

print("Calculating depth with the highest number of sites to truncate x axis...")
max_depth <- depth_global %>%
  filter(count_good_global == max(count_good_global)) %>%
  pull(depth)

print("Plotting genome-wide depth distribution...")
plot_global <- ggplot()+
  geom_line(data = depth_global, aes(x=depth,y = count_raw_global, color = "All alignments"))+
  geom_line(data = depth_global, aes(x=depth,y = count_good_global, color = "Good quality alignments"))+ 
  scale_y_log10(name = "Number of Sites", labels = comma)+
  scale_x_continuous(name = "Depth (X)", labels = comma, n.breaks = 10)+
  labs(subtitle = "Log10 scale in Y-axis")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        legend.position = "none")

plot_truncated_log <- ggplot()+
  geom_line(data = depth_global, aes(x=depth,y = count_raw_global, color = "All alignments"))+ 
  geom_line(data = depth_global, aes(x=depth,y = count_good_global, color = "Good quality alignments"))+ 
  scale_y_log10(name = "Number of Sites", labels = comma)+
  scale_x_continuous(name = "Depth (X)", labels = comma, n.breaks = 10, limits = c(0,max_depth*10))+
  labs(subtitle = "Log10 scale in Y-axis, truncated X-axis") +
  theme_bw()+
  theme(axis.title.x = element_blank(),
        legend.position = "none")

plot_truncated <- ggplot()+
  geom_line(data = depth_global, aes(x=depth,y = count_raw_global, color = "All alignments"))+ 
  geom_line(data = depth_global, aes(x=depth,y = count_good_global, color = "Good quality alignments"))+ 
  scale_y_continuous(name = "Number of Sites", labels = comma)+
  scale_x_continuous(name = "Depth (X)", labels = comma, n.breaks = 10, limits = c(0,max_depth*10))+
  labs(subtitle = "Truncated X-axis")+
  theme_bw() +
  theme(legend.position = "bottom")+
  guides(color = guide_legend(nrow = 1, title = "Quality"))


combined <- plot_global / plot_truncated_log / plot_truncated 
combined <- combined +
  plot_annotation(title = "Depth Distribution of Whole Genome by Quality of Read Alignments",
                  subtitle = paste("Lineage:", lineage_name, " Sample:", sample,"Strain:", strain_name, sep = " ")) &
    theme(plot.title = element_text(hjust = 0.5))

print("Plotting depth distribution by chromosome...")
by_chrom <- ggplot(depth)+
  geom_line(aes(x=depth, y = count_good, color = chromosome))+ 
  scale_y_log10(name = "Number of Sites", labels = comma)+
  scale_x_continuous(name = "Depth (X)", labels = comma, n.breaks = 10)+
  labs(subtitle = "Log10 scale in Y-axis")+
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

by_chrom_log <- ggplot(depth)+
  geom_line(aes(x=depth, y = count_good, color = chromosome))+ 
  scale_y_log10(name = "Number of Sites", labels = comma)+
  scale_x_continuous(name = "Depth (X)", labels = comma, n.breaks = 10, limits = c(0,max_depth*10))+
  labs(subtitle = "Log10 scale in Y-axis, truncated X-axis")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        legend.position = "none")

by_chrom_truncated <- ggplot(depth)+
  geom_line(aes(x=depth, y = count_good, color = chromosome))+ 
  scale_y_continuous(name = "Number of Sites", labels = comma)+
  scale_x_continuous(name = "Depth (X)", labels = comma, n.breaks = 10, limits = c(0,max_depth*10))+
  labs(subtitle = "Truncated X-axis")+
  theme_bw() +
  theme(legend.position = "bottom", legend.direction= "horizontal") +
  guides(color = guide_legend(nrow = 1, title = "Chromosome"))


plot_chrom <- by_chrom / by_chrom_log / by_chrom_truncated
plot_chrom <- plot_chrom +
  plot_annotation(title = "Depth Distribution of Good Quality Alignments of each Chromosome",
                  subtitle = paste("Lineage:", lineage_name, " Sample:", sample, "Strain:", strain_name, sep = " ")) &
    theme(plot.title = element_text(hjust = 0.5))

print("Saving plots...")
ggsave(output_chrom, plot = plot_chrom, units = "in", height = gheight, width = gwidth, dpi = gdpi, scale = gscale)
ggsave(output_global, plot = combined,  units = "in", height = gheight, width = gwidth, dpi = gdpi, scale = gscale)

print("Done!")







