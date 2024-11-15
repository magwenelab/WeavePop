log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggnewscale))

print("Reading files...")
raw<- read.delim(snakemake@input[[1]], header = FALSE, col.names = c("accession", "start", "end", "mapq"), stringsAsFactors = TRUE)
cnv <- read.delim(snakemake@input[[2]], sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
repeats_table <- read.delim(snakemake@input[[3]], sep= "\t", header = FALSE, col.names = c("accession", "start", "end", "repeat_type"), stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
chrom_names <- read.csv(snakemake@input[[4]], header = FALSE, col.names = c("lineage", "accession", "chromosome"))
loci <- read.delim(snakemake@input[[5]], header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
sample <- snakemake@wildcards$sample

print("Filtering chromosome names...")
chrom_names <- chrom_names %>%
  filter(accession %in% unique(raw$accession) )

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
raw <- left_join(raw, chrom_names, by = "accession")
loci <- loci %>% 
  select(accession, start, loci)%>%
  filter(accession %in% unique(raw$accession))
lineage <- levels(as.factor(raw$lineage))

cnv <- left_join(cnv, chrom_names, by = "accession")
cnv$cnv <- str_to_title(cnv$cnv)
cnv$cnv <- as.factor(cnv$cnv)

feature <- cnv %>%
  select(accession_chromosome, chromosome, start, end, cnv)

repeats<- left_join(repeats_table, chrom_names, by = "accession")%>%
  select(accession_chromosome, chromosome, start, end, repeat_type)
repeats$repeat_type <- ifelse(repeats$repeat_type == "Simple_repeat", "Simple repeat", "Others")
repeats$repeat_type <- factor(repeats$repeat_type, levels = c("Simple repeat", "Others"))

print("Getting plot parameters...")
s_lim <- max(raw$mapq) + 40
set2 <- rev(brewer.pal(8, "Set2")[1:6])
s_colors <- set2[1:nlevels(feature$cnv)]
r_lim <- max(raw$mapq) + 60
r_colors <- colorRampPalette(brewer.pal(12, "Paired"))(nlevels(repeats$repeat_type))
raw_color <- "black"
palette1 <- brewer.pal(n = 7, name = "Dark2")
palette2 <- brewer.pal(n = 7, name = "Set2")
combined_palette <- c(palette1, palette2)

print("Creating labeller function...")
my_labeller <- function(value){
  new_value <- sapply(strsplit(as.character(value), "xxx"), head, 1)
  return(new_value)
}
print("Plotting chromosome MAPQ...")
plot <- ggplot()+
  coord_cartesian(ylim=c(0,r_lim + 20), xlim= c(0,max(raw$end)))+
  geom_point(data = raw, aes(x= start, y = mapq), size = 0.5 , color = raw_color)+
  geom_segment(data = repeats, aes(x = start, xend = end, y = r_lim, yend = r_lim, color = repeat_type), linewidth = 2)+
    scale_color_manual(name = "Type of repetitive sequence", values = r_colors)+
    guides(color = guide_legend(order=1))+
    new_scale_color()+
  geom_segment(data = feature, aes(x = start, xend = end, y = s_lim, yend = s_lim, color = cnv), linewidth = 2)+
    scale_color_manual(name = "Copy-number variant", values = s_colors)+
    guides(color = guide_legend(order=2))+
    new_scale_color()+
  facet_wrap(~accession_chromosome, strip.position = "right", ncol = 2, labeller = as_labeller(my_labeller)) +
  scale_x_continuous(name = "Position (bp) ", labels = comma)+
  labs(title = paste("Lineage: ", lineage," Sample: ", sample,  sep = " "), y = "Mapping quality (phred score)")+
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 2))+
  scale_y_continuous(breaks = c(0, 20, 40, 60))

print("Adding loci to plot if available...")
y_loci = max(raw$mapq) + 20
if (nrow(loci) != 0){
  loci <- loci %>%
    filter(lineage %in% lineage)
  loci_chrom <- left_join(loci, chrom_names, by = "accession")
  loci_chrom <- loci_chrom %>% filter(lineage %in% lineage)
  loci_colors <- combined_palette[1:nlevels(loci_chrom$loci)]
  plot <- plot +
    geom_point(data = loci_chrom, aes(x= start, y = y_loci, color = loci), size = 2, shape = 19)+
    scale_color_manual(name = "Features", values = loci_colors)
}

print("Saving plot...")
ggsave(snakemake@output[[1]], plot = plot, units = "in", height = 9, width = 16, dpi = 600)

print("Done!")
