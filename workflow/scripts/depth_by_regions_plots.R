log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(RColorBrewer))

print("Reading files")
coverage_regions <- read.delim(snakemake@input[[1]], sep= "\t", col.names = c("Accession", "Start", "End", "Depth", "Norm_Depth", "Smooth_Depth"), stringsAsFactors = TRUE, na = c("", "N/A"))
struc_vars <- read.delim(snakemake@input[[2]], sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
repeats_table <- read.delim(snakemake@input[[3]], sep= "\t", header = FALSE, col.names = c("Accession", "Start", "End", "Repeat_type"), stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
loci_table <- read.delim(snakemake@input[[4]], header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
chrom_names <- read.csv(snakemake@input[[5]], sep = ",", header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"), stringsAsFactors = TRUE, na = c("", "N/A"))
sample <- snakemake@wildcards$sample

# coverage_regions <- read.delim("results/samples/mosdepth/PMY3493/good_regions_coverage.tsv", sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
# struc_vars <- read.delim("results/samples/mosdepth/PMY3493/good_structural_variants.tsv", sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
# repeats_table <- read.delim("results/references/VNI_VNIV/repeats/VNI_VNIV_repeats.bed", sep= "\t", header = FALSE, col.names = c("Accession", "Start", "End", "Repeat_type"), stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
# loci_table <- read.delim("results/dataset/files/loci_to_plot.tsv", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
# chrom_names <- read.csv("config/chromosome_names.csv", header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"), stringsAsFactors = TRUE, na = c("", "N/A"))
# sample <- "PMY3493"
print("Get chromosome names")
chrom_names <- chrom_names %>%
  filter(Accession %in% unique(struc_vars$Accession) )
chrom_names['Accession_Chromosome'] <- paste(chrom_names$Chromosome, chrom_names$Accession, sep = "xxx")
unique_levels <- unique(chrom_names$Accession_Chromosome)
new_order <- c(rbind(matrix(unique_levels, nrow = 2, byrow = TRUE)))
chrom_names$Accession_Chromosome <- factor(chrom_names$Accession_Chromosome, levels = new_order)

print("Rearrange loci data")
loci_sample <- loci_table %>% 
    select(Accession = seq_id, Start = start, End = end,Loci)%>%
    filter(Accession %in% coverage_regions$Accession)%>%
    mutate(Track = "Loci")%>%
    droplevels()
loci <- left_join(loci_sample, chrom_names, by = c("Accession"))%>%
  select(Accession_Chromosome, Chromosome, Track, Start, End, Loci)
dark2 <- brewer.pal(8, "Dark2")[1:6]
l_colors <- dark2[1:nlevels(loci$Loci)]

coverage_regions <- left_join(coverage_regions, chrom_names, by = "Accession")
coverage <- coverage_regions %>%
  select(Accession_Chromosome, Chromosome, Start, End, Coverage = Norm_Depth)%>%
  mutate(Track = "Coverage", .after = Chromosome)
topCov <- quantile(coverage$Coverage, 0.75) * 3
coverage$Coverage<- ifelse(coverage$Coverage >= topCov, topCov, coverage$Coverage)
l_lim <- topCov 

smooth <- coverage_regions %>%
  select(Accession_Chromosome, Chromosome, Start, End, Smooth = Smooth_Depth)%>%
  mutate(Track = "Smooth", .after = Chromosome)


struc_vars <- left_join(struc_vars, chrom_names, by = "Accession")
structure <- struc_vars %>%
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

lineage <- unique(coverage_regions$Lineage)

# variants <- read.delim(snakemake@input[[4]], sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
# variants <- read.delim("results/dataset/snps/VNI_variants.tsv", sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A", "NA"))

# variants <- variants %>% 
#   select(VAR, Accession = CHROM, Start = POS, Sample = all_of(sample), Gene = any_of(c("ID", "gene_id", "locus_tag")), Description = any_of(c("description", "product"))) %>%
#   mutate(Track = "Variants")%>%
#   filter(Sample == 1)%>%
#   droplevels()
# variants <- left_join(variants, chrom_names, by = "Accession")
# v_lim <- topCov + 3

my_labeller <- function(value){
  new_value <- sapply(strsplit(as.character(value), "xxx"), head, 1)
  return(new_value)
}

# Coverage plot
c <- ggplot()+
  coord_cartesian(ylim= c(0,r_lim +1), xlim = c(0, max(coverage$End)))+
  geom_hline(yintercept = 1, color = "darkgray", linetype = 2)+
  geom_hline(yintercept = 2, color = "darkgray", linetype = 2)+
  geom_col(data = coverage, aes(x=Start, y = Coverage), color = "black")+
    scale_x_continuous(name = "Position (bp) ", labels = comma)+
  geom_segment(data = repeats, aes(x = Start, xend = End, y = r_lim, yend = r_lim, color = Repeat_type), linewidth = 2)+
    scale_color_manual(name = "Type of repetitive sequence", values = r_colors)+
    guides(color = guide_legend(order=1))+
    new_scale_color()+
  geom_segment(data = structure, aes(x = Start, xend = End, y = s_lim, yend = s_lim, color = Structure), linewidth = 2)+
    scale_color_manual(name = "Copy number variants", values = s_colors)+
    guides(color = guide_legend(order=2))+
    new_scale_color()+
  geom_point(data = loci, aes(x=Start, y = l_lim, color = Loci))+  
    scale_color_manual(name = "Features", values = l_colors)+
    guides(color = guide_legend(order=3))+
  # geom_point(data = variants, aes(x=Start, y = v_lim),shape = 24 , size = 2)+
  facet_wrap(~Accession_Chromosome, strip.position = "right", ncol = 2, labeller = as_labeller(my_labeller)) +
  labs(y = "Normalized coverage", title = paste("Lineage:", lineage, " Sample:", sample, sep = " "))+
  scale_y_continuous(breaks = c(1, 2)) +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 2))

ggsave(snakemake@output[[1]], c, height =9, width = 16, dpi = 600)

