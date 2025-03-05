log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scales))

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

print("Joining and arranging data...")
depth <- left_join(depth, chrom_names, by = "accession")
lineage <- levels(as.factor(depth$lineage))

print("Getting plot parameters...")
raw_color = "gray50"
good_color = "black" 
color_quality = c("Good quality alignments" = good_color, "All alignments" = raw_color)

print("Plotting depth distribution by chromosome...")
plot_chrom <- ggplot(depth, aes(x=depth))+
  geom_col(aes(y = count_raw, fill= "All alignments"))+ 
  geom_col(aes(y = count_good, fill= "Good quality alignments"))+ 
  facet_wrap(~chromosome,ncol = 2)+
  scale_y_log10(name = "Number of Sites", labels = comma)+
  scale_x_continuous(name = "Depth (X) ", labels = comma, n.breaks = 10)+
  scale_fill_manual(name= "Alignment quality", values= color_quality)+
  # theme(legend.position="none")+
  labs(title = paste("Lineage:",lineage, " Sample:", sample,  sep = " "))+
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major = element_line(color = "gray95"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 2))



print("Calculating global depth distribution...")
depth_global <- depth %>%
  select(depth, count_good, count_raw)%>%
  group_by(depth)%>%
  summarize(count_good_global = sum(count_good), count_raw_global = sum(count_raw))%>%
  ungroup()

print("Plotting genome-wide depth distribution...")
plot_global <- ggplot()+
  geom_col(data = depth_global, aes(x=depth, y = count_raw_global, fill= "All alignments"))+ 
  geom_col(data = depth_global, aes(x=depth,y = count_good_global, fill= "Good quality alignments"))+ 
  scale_y_log10(name = "Number of Sites", labels = comma)+
  scale_x_continuous(name = "Depth (X) ", labels = comma, n.breaks = 10)+
  scale_fill_manual(name= "Alignment quality", values= color_quality)+
  theme_bw()+
  theme(legend.position="none",
    axis.title = element_blank(),
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
  # labs(title = paste("Lineage:",lineage, " Sample:", sample,  sep = " "))

# getting mode
max_depth <- depth_global %>%
  filter(count_good_global == max(count_good_global)) %>%
  pull(depth)

# Truncated plot
plot_truncated <- ggplot()+
  geom_col(data = depth_global, aes(x=depth, y = count_raw_global, fill= "All alignments"))+ 
  geom_col(data = depth_global, aes(x=depth,y = count_good_global, fill= "Good quality alignments"))+ 
  scale_y_log10(name = "Number of Sites", labels = comma)+
  scale_x_continuous(name = "Depth (X) ", labels = comma, n.breaks = 10, limits = c(0,max_depth*10))+
  scale_fill_manual(name= "Alignment quality", values= color_quality)+
  theme(legend.position="none")+
  labs(title = paste("Lineage:",lineage, " Sample:", sample,  sep = " "))+
  theme_bw()

overlay_grob <- ggplotGrob(plot_global)

# Add the overlay plot to the main plot
combined_plot <- plot_truncated +
  annotation_custom(grob = overlay_grob,
      xmin = 0.5 * max_depth * 10,
      xmax = max_depth * 10, 
      ymin = 0.5 * max(log10(depth_global$count_good_global)),
      ymax = max(log10(depth_global$count_good_global)))

print("Saving plots...")
ggsave(snakemake@output[[1]], plot = plot_chrom, units = "in", height = 4.5, width = 8, dpi = 600)
# ggsave(snakemake@output[[2]], plot = plot_global,  units = "in", height = 4.5, width = 8, dpi = 600)
ggsave(snakemake@output[[2]], plot = combined_plot,  units = "in", height = 4.5, width = 8, dpi = 600)


print("Done!")




