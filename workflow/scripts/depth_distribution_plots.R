log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(patchwork))


print("Reading files and joining data with chromosome names...")
sample <- snakemake@wildcards$sample
depth<- read.table(snakemake@input[[1]], header = TRUE, stringsAsFactors = TRUE, sep = "\t")
chrom_names <- read.csv(snakemake@input[[2]], header = TRUE, col.names = c("lineage", "accession", "chromosome"))

# sample <- "ERS542301"
# depth <- read.table("/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Ashton/results/04.Intermediate_files/01.Samples/depth_quality/ERS542301/depth_distribution.tsv", header = TRUE, stringsAsFactors = TRUE, sep = "\t")
# chrom_names <- read.csv("/FastData/czirion/Crypto_Diversity_Pipeline/Crypto_Ashton/results/02.Dataset/chromosomes.csv", header = TRUE, col.names = c("lineage", "accession", "chromosome"))

print("Adding 0 to null depth values...")
depth[is.na(depth)] <- 0

print("Filtering chromosome names...")
chrom_names <- chrom_names %>%
  filter(accession %in% unique(depth$accession) )

print("Ordering chromosome names...")
chrom_names['accession_chromosome'] <- paste(chrom_names$chromosome, chrom_names$accession, sep = "xxx")
unique_levels <- unique(chrom_names$accession_chromosome)
new_order <- c(rbind(matrix(unique_levels, nrow = 2, byrow = TRUE)))
# If unique levels is an odd number, the last element is removed
if (length(unique_levels) %% 2 != 0){
  new_order <- new_order[-length(new_order)]
}
chrom_names$accession_chromosome <- factor(chrom_names$accession_chromosome, levels = new_order)
chrom_names$chromosome <- factor(chrom_names$chromosome, levels = unique(chrom_names$chromosome))

print("Joining and arranging data...")
depth <- left_join(depth, chrom_names, by = "accession")
lineage <- levels(as.factor(depth$lineage))

print("Getting plot parameters...")
raw_color = "gray50"
good_color = "black" 
color_quality = c("Good quality alignments" = good_color, "All alignments" = raw_color)

print("Calculating global depth distribution...")
depth_global <- depth %>%
  select(depth, count_good, count_raw)%>%
  group_by(depth)%>%
  summarize(count_good_global = sum(count_good), count_raw_global = sum(count_raw))%>%
  ungroup()

print("Calculating depth witht the highest number of sites to truncate x axis...")
max_depth <- depth_global %>%
  filter(count_good_global == max(count_good_global)) %>%
  pull(depth)

print("Plotting genome-wide depth distribution...")
plot_global <- ggplot()+
  geom_line(data = depth_global, aes(x=depth,y = count_good_global), fill = good_color)+ 
  scale_y_log10(name = "Number of Sites", labels = comma)+
  scale_x_continuous(name = "Depth (X)", labels = comma, n.breaks = 10)+
  theme_bw()+
  theme(axis.title.x = element_blank())+
  labs(title = "Good quality alignments", subtitle = "Log10 scale in Y-axis")

plot_truncated_log <- ggplot()+
  geom_line(data = depth_global, aes(x=depth,y = count_good_global), fill = good_color)+ 
  scale_y_log10(name = "Number of Sites", labels = comma)+
  scale_x_continuous(name = "Depth (X)", labels = comma, n.breaks = 10, limits = c(0,max_depth*10))+
  theme_bw()+
  theme(axis.title.x = element_blank())+
  labs(subtitle = "Log10 scale in Y-axis, truncated X-axis")

plot_truncated <- ggplot()+
  geom_line(data = depth_global, aes(x=depth,y = count_good_global), fill = good_color)+ 
  scale_y_continuous(name = "Number of Sites", labels = comma)+
  scale_x_continuous(name = "Depth (X)", labels = comma, n.breaks = 10, limits = c(0,max_depth*10))+
  theme_bw()+
  labs(subtitle = "Truncated X-axis")

plot_global_raw <- ggplot()+
  geom_line(data = depth_global, aes(x=depth,y = count_raw_global), fill = raw_color)+ 
  scale_y_log10(name = "Number of Sites", labels = comma)+
  scale_x_continuous(name = "Depth (X)", labels = comma, n.breaks = 10)+
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
  labs(title = "All alignments", subtitle = "Log10 scale in Y-axis")

plot_truncated_log_raw <- ggplot()+
  geom_line(data = depth_global, aes(x=depth,y = count_raw_global), fill = raw_color)+ 
  scale_y_log10(name = "Number of Sites", labels = comma)+
  scale_x_continuous(name = "Depth (X)", labels = comma, n.breaks = 10, limits = c(0,max_depth*10))+
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
  labs(subtitle = "Log10 scale in Y-axis, truncated X-axis")

plot_truncated_raw <- ggplot()+
  geom_line(data = depth_global, aes(x=depth,y = count_raw_global), fill = raw_color)+ 
  scale_y_continuous(name = "Number of Sites", labels = comma)+
  scale_x_continuous(name = "Depth (X)", labels = comma, n.breaks = 10, limits = c(0,max_depth*10))+
  theme_bw()+
  theme(axis.title.y = element_blank())+
  labs(subtitle = "Truncated X-axis")

combined <- plot_global / plot_truncated_log / plot_truncated | plot_global_raw / plot_truncated_log_raw / plot_truncated_raw
combined <- combined +
  plot_annotation(title = paste("Depth Distribution for Sample", sample, " of Lineage", lineage)) &
    theme(plot.title = element_text(hjust = 0.5))


print("Plotting depth distribution by chromosome...")
by_chrom <- ggplot(depth)+
  geom_line(aes(x=depth, y = count_good, color = chromosome))+ 
  scale_y_continuous(name = "Number of Sites", labels = comma)+
  scale_x_continuous(name = "Depth (X)", labels = comma, n.breaks = 10, limits = c(0,max_depth*10))+
  labs(subtitle = "Truncated X-axis")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        legend.position = "none")

by_chrom_log <- ggplot(depth)+
  geom_line(aes(x=depth, y = count_good, color = chromosome))+ 
  scale_y_log10(name = "Number of Sites", labels = comma)+
  scale_x_continuous(name = "Depth (X)", labels = comma, n.breaks = 10, limits = c(0,max_depth*10))+
  labs(subtitle = "Log10 scale in Y-axis, truncated X-axis")+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction= "horizontal") +
  guides(color = guide_legend(override.aes = list(size = 5), nrow = 1))

plot_chrom <- by_chrom / by_chrom_log
plot_chrom <- plot_chrom +
  plot_annotation(title = paste("Depth Distribution for Sample", sample, " of Lineage", lineage)) &
    theme(plot.title = element_text(hjust = 0.5))

print("Saving plots...")
ggsave(snakemake@output[[1]], plot = plot_chrom, units = "in", height = 4.5, width = 8, dpi = 600, scale = 1.5)
ggsave(snakemake@output[[2]], plot = combined,  units = "in", height = 4.5, width = 8, dpi = 600, scale = 2)


print("Done!")







