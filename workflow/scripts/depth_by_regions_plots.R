log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(RColorBrewer))

print("Reading files")
depth_regions <- read.delim(snakemake@input[[1]], sep= "\t", col.names = c("Accession", "Start", "End", "Depth", "Norm_Depth", "Smooth_Depth"), stringsAsFactors = TRUE, na = c("", "N/A"))
cnv <- read.delim(snakemake@input[[2]], sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
repeats_table <- read.delim(snakemake@input[[3]], sep= "\t", header = FALSE, col.names = c("Accession", "Start", "End", "Repeat_type"), stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
loci_table <- read.delim(snakemake@input[[4]], header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
chrom_names <- read.csv(snakemake@input[[5]], sep = ",", header = TRUE, col.names = c("Lineage", "Accession", "Chromosome"), stringsAsFactors = TRUE, na = c("", "N/A"))
sample <- snakemake@wildcards$sample

# depth_regions <- read.delim("results_230724/samples/mosdepth/PMY3315/depth_by_regions.tsv", sep= "\t", col.names = c("Accession", "Start", "End", "Depth", "Norm_Depth", "Smooth_Depth"), stringsAsFactors = TRUE, na = c("", "N/A"))
# cnv <- read.delim("results_230724/samples/cnv/PMY3315/cnv_calls.tsv", sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
# repeats_table <- read.delim("results_230724/references/Cdeneoformans/repeats/Cdeneoformans_repeats.bed", sep= "\t", header = FALSE, col.names = c("Accession", "Start", "End", "Repeat_type"), stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
# loci_table <- read.delim("results_230724/dataset/files/loci_to_plot.tsv", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
# chrom_names <- read.csv("config/chromosome_names.csv", sep = ",", header = TRUE, col.names = c("Lineage", "Accession", "Chromosome"), stringsAsFactors = TRUE, na = c("", "N/A"))
# sample <- "PMY3315"
print("Getting chromosome names")
chrom_names <- chrom_names %>%
  filter(Accession %in% unique(depth_regions$Accession) )
chrom_names['Accession_Chromosome'] <- paste(chrom_names$Chromosome, chrom_names$Accession, sep = "xxx")
unique_levels <- unique(chrom_names$Accession_Chromosome)
new_order <- c(rbind(matrix(unique_levels, nrow = 2, byrow = TRUE)))
# If unique levels is an odd number, the last element is removed
if (length(unique_levels) %% 2 != 0){
  new_order <- new_order[-length(new_order)]
}
chrom_names$Accession_Chromosome <- factor(chrom_names$Accession_Chromosome, levels = new_order)

print("Joining and arranging data")
depth_regions <- left_join(depth_regions, chrom_names, by = "Accession")
depth <- depth_regions %>%
  select(Accession_Chromosome, Chromosome, Start, End, Depth = Norm_Depth)%>%
  mutate(Track = "Depth", .after = Chromosome)
topCov <- quantile(depth$Depth, 0.75) * 3
depth$Depth<- ifelse(depth$Depth >= topCov, topCov, depth$Depth)

smooth <- depth_regions %>%
  select(Accession_Chromosome, Chromosome, Start, End, Smooth = Smooth_Depth)%>%
  mutate(Track = "Smooth", .after = Chromosome)

cnv$Structure <- str_to_title(cnv$Structure)
cnv$Structure <- as.factor(cnv$Structure)
cnv <- left_join(cnv, chrom_names, by = "Accession")
structure <- cnv %>%
  select(Accession_Chromosome, Chromosome, Start, End, Structure)%>%
  mutate(Track = "Copy_Number_Variants")
s_lim <- topCov + 1
set2 <- rev(brewer.pal(8, "Set2")[1:6])
s_colors <- set2[1:nlevels(structure$Structure)]

repeats<- left_join(repeats_table, chrom_names, by = "Accession")%>%
  select(Accession_Chromosome, Chromosome, Start, End, Repeat_type)%>%
  mutate(Track = "Repeats", .after= Chromosome)

repeats$Repeat_type <- ifelse(repeats$Repeat_type == "Simple_repeat", "Simple repeat", "Others")
repeats$Repeat_type <- factor(repeats$Repeat_type, levels = c("Simple repeat", "Others"))
r_lim <- topCov + 2
r_colors <- colorRampPalette(brewer.pal(12, "Paired"))(nlevels(repeats$Repeat_type))

lineage <- unique(depth_regions$Lineage)

# variants <- read.delim(snakemake@input[[4]], sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
# variants <- read.delim("results/dataset/snps/VNI_variants.tsv", sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A", "NA"))

# variants <- variants %>% 
#   select(VAR, Accession = CHROM, Start = POS, Sample = all_of(sample), Gene = any_of(c("ID", "gene_id", "locus_tag")), Description = any_of(c("description", "product"))) %>%
#   mutate(Track = "Variants")%>%
#   filter(Sample == 1)%>%
#   droplevels()
# variants <- left_join(variants, chrom_names, by = "Accession")
# v_lim <- topCov + 3

print("Making labeler function")
my_labeller <- function(value){
  new_value <- sapply(strsplit(as.character(value), "xxx"), head, 1)
  return(new_value)
}

print("Plotting")
# Depth plot
c <- ggplot()+
  coord_cartesian(ylim= c(0,r_lim +1), xlim = c(0, max(depth$End)))+
  geom_hline(yintercept = 1, color = "darkgray", linetype = 2)+
  geom_hline(yintercept = 2, color = "darkgray", linetype = 2)+
  geom_col(data = depth, aes(x=Start, y = Depth), color = "black")+
    scale_x_continuous(name = "Position (bp) ", labels = comma)+
  geom_segment(data = repeats, aes(x = Start, xend = End, y = r_lim, yend = r_lim, color = Repeat_type), linewidth = 2)+
    scale_color_manual(name = "Type of repetitive sequence", values = r_colors)+
    guides(color = guide_legend(order=1))+
    new_scale_color()+
  geom_segment(data = structure, aes(x = Start, xend = End, y = s_lim, yend = s_lim, color = Structure), linewidth = 2)+
    scale_color_manual(name = "Copy number variants", values = s_colors)+
    guides(color = guide_legend(order=2))+
    new_scale_color()+
  # geom_point(data = variants, aes(x=Start, y = v_lim),shape = 24 , size = 2)+
  facet_wrap(~Accession_Chromosome, strip.position = "right", ncol = 2, labeller = as_labeller(my_labeller)) +
  labs(y = "Normalized depth", title = paste("Lineage:", lineage, " Sample:", sample, sep = " "))+
  scale_y_continuous(breaks = c(1, 2)) +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 2))

print("Adding loci data if available")
if (nrow(loci_table)!= 0){
print("Rearrange loci data")
loci_sample <- loci_table %>% 
    select(Accession = seq_id, Start = start, End = end,Loci)%>%
    filter(Accession %in% depth_regions$Accession)%>%
    mutate(Track = "Loci")%>%
    droplevels()
loci <- left_join(loci_sample, chrom_names, by = c("Accession"))%>%
  select(Accession_Chromosome, Chromosome, Track, Start, End, Loci)
dark2 <- brewer.pal(8, "Dark2")[1:6]
l_colors <- dark2[1:nlevels(loci$Loci)]
l_lim <- topCov 

c <- c +  geom_point(data = loci, aes(x=Start, y = l_lim, color = Loci))+  
    scale_color_manual(name = "Features", values = l_colors)+
    guides(color = guide_legend(order=3))
}

print("Saving plot")
ggsave(snakemake@output[[1]], c, height =9, width = 16, dpi = 600)
print("Done!")
