log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(RColorBrewer))

print("Reading files...")
depth_windows <- read.delim(snakemake@input[[1]], sep= "\t", col.names = c("accession", "start", "end", "depth", "norm_depth", "smooth_depth"), stringsAsFactors = TRUE, na = c("", "N/A"))
cnv <- read.delim(snakemake@input[[2]], sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
repeats_table <- read.delim(snakemake@input[[3]], sep= "\t", header = FALSE, col.names = c("accession", "start", "end", "repeat_type"), stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
loci_table <- read.delim(snakemake@input[[4]], header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
chrom_names <- read.csv(snakemake@input[[5]], sep = ",", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
sample <- snakemake@wildcards$sample
metadata <- read.delim(snakemake@input[[6]], sep = ",", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))

print("Obtaining lineage of sample...")

lineage_name <- as.character(metadata$lineage[metadata$sample == sample])

print("Filtering chromosome names...")
chrom_names <- chrom_names %>%
  filter(lineage == lineage_name)

print("Ordering chromosome names...")
chrom_names['accession_chromosome'] <- paste(chrom_names$chromosome, chrom_names$accession, sep = "xxx")
unique_levels <- unique(chrom_names$accession_chromosome)
new_order <- c(rbind(matrix(unique_levels, nrow = 2, byrow = TRUE)))
# If unique levels is an odd number, the last element is removed
if (length(unique_levels) %% 2 != 0){
  new_order <- new_order[-length(new_order)]
}
chrom_names$accession_chromosome <- factor(chrom_names$accession_chromosome, levels = new_order)

print("Arranging depth data...")
depth_windows <- left_join(depth_windows, chrom_names, by = "accession")
depth <- depth_windows %>%
  select(accession_chromosome, chromosome, start, end, depth = norm_depth)

print("Truncating extremele high depth values...")
topcov <- quantile(depth$depth, 0.75) * 3
depth$depth<- ifelse(depth$depth >= topcov, topcov, depth$depth)

print("Arranging CNV data...")
cnv$cnv <- str_to_title(cnv$cnv)
cnv$cnv <- as.factor(cnv$cnv)
cnv <- left_join(cnv, chrom_names, by = "accession")

feature <- cnv %>%
  select(accession_chromosome, chromosome, start, end, cnv)

print("Arranging repeats data...")
repeats<- left_join(repeats_table, chrom_names, by = "accession")%>%
  select(accession_chromosome, chromosome, start, end, repeat_type)

repeats$repeat_type <- ifelse(repeats$repeat_type == "Simple_repeat", "Simple repeat", "Others")
repeats$repeat_type <- factor(repeats$repeat_type, levels = c("Simple repeat", "Others"))

print("Getting plot parameters...")
s_lim <- topcov + 1
set2 <- rev(brewer.pal(8, "Set2")[1:6])
s_colors <- set2[1:nlevels(feature$cnv)]
r_lim <- topcov + 2
r_colors <- colorRampPalette(brewer.pal(12, "Paired"))(nlevels(repeats$repeat_type))

print("Making labeler function...")
my_labeller <- function(value){
  new_value <- sapply(strsplit(as.character(value), "xxx"), head, 1)
  return(new_value)
}

print("Plotting depth by windows...")
c <- ggplot()+
  coord_cartesian(ylim= c(0,r_lim +1), xlim = c(0, max(depth$end)))+
  geom_hline(yintercept = 1, color = "darkgray", linetype = 2)+
  geom_hline(yintercept = 2, color = "darkgray", linetype = 2)+
  geom_col(data = depth, aes(x=start, y = depth), color = "black")+
    scale_x_continuous(name = "Position (bp) ", labels = comma)+
  geom_segment(data = repeats, aes(x = start, xend = end, y = r_lim, yend = r_lim, color = repeat_type), linewidth = 2)+
    scale_color_manual(name = "Type of Repetitive\nSequences", values = r_colors)+
    guides(color = guide_legend(order=1))+
    new_scale_color()+
  geom_segment(data = feature, aes(x = start, xend = end, y = s_lim, yend = s_lim, color = cnv), linewidth = 2)+
    scale_color_manual(name = "Copy-Number\nVariants", values = s_colors)+
    guides(color = guide_legend(order=2))+
    new_scale_color()+
  facet_wrap(~accession_chromosome, strip.position = "right", ncol = 2, labeller = as_labeller(my_labeller)) +
  labs(y = "Normalized Depth",
      title = "Normalized Depth of Windows Along Chromosomes", 
      subtitle = paste("Lineage:", lineage_name, " Sample:", sample, sep = " "))+
  scale_y_continuous(breaks = c(1, 2)) +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 2))

print("Adding loci data if available...")
if (nrow(loci_table)!= 0){
  loci_sample <- loci_table %>% 
    filter(lineage == lineage_name) %>%
    select(accession, start , end , loci)%>%
    droplevels()
  loci <- left_join(loci_sample, chrom_names, by = c("accession"))%>%
    select(accession_chromosome, chromosome, start, end, loci)
  dark2 <- brewer.pal(8, "Dark2")[1:6]
  l_colors <- dark2[1:nlevels(loci$loci)]
  l_lim <- topcov 
  print("Adding loci to plot...")
  c <- c +  
      geom_point(data = loci, aes(x=start, y = l_lim, color = loci))+  
      scale_color_manual(name = "Features", values = l_colors)+
      guides(color = guide_legend(order=3))
}

print("Saving plot...")
ggsave(snakemake@output[[1]], c, height =9, width = 16, dpi = 600)

print("Done!")
