suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(RColorBrewer))


sample <- "results/samples/mosdepth/SRS8318899/smooth_coverage_regions.tsv"
Split <- str_split(sample, "/")
sample <- Split[[1]][length(Split[[1]])-1]

good_stats_regions <- read.delim("results/samples/mosdepth/SRS8318899/smooth_good_stats_regions.tsv", sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))
struc_vars <- read.delim("results/samples/mosdepth/SRS8318899/ploidy_table.tsv", sep= "\t", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
repeats_table <- read.delim("results/references/VNI/repeats/05_full/VNI.bed", sep= "\t", header = FALSE, col.names = c("Accession", "Start", "End", "Repeat_type"), stringsAsFactors = TRUE, na = c("", "N/A", "NA"))
loci_table <- read.delim("results/dataset/loci_to_plot.tsv", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))

chrom_names <- good_stats_regions %>%
    select(Accession, Chromosome, Lineage)%>%
    distinct()
lineage <- levels(as.factor(chrom_names$Lineage))

loci_sample <- loci_table %>% 
    select(Accession = seq_id, Start = start, End = end,Loci)%>%
    filter(Accession %in% good_stats_regions$Accession)%>%
    mutate(Track = "Loci")%>%
    droplevels()
loci <- left_join(loci_sample, chrom_names, by = c("Accession"))%>%
  select(Chromosome, Track, Start, End, Loci)
dark2 <- brewer.pal(8, "Dark2")[1:6]
l_colors <- dark2[1:nlevels(loci$Loci)]


coverage <- good_stats_regions %>%
  select(Chromosome, Start, End, Coverage = Norm_Median)%>%
  mutate(Track = "Coverage", .after = Chromosome)
topCov <- quantile(coverage$Coverage, 0.75) * 3
coverage$Coverage<- ifelse(coverage$Coverage >= topCov, topCov, coverage$Coverage)
l_lim <- topCov 

structure <- struc_vars %>%
  select(Chromosome, Start, End, Structure)%>%
  mutate(Track = "Structural_Variants")
s_lim <- topCov + 1
set2 <- rev(brewer.pal(8, "Set2")[1:6])
s_colors <- set2[1:nlevels(structure$Structure)]

repeats<- left_join(repeats_table, chrom_names, by = "Accession")%>%
  select(Chromosome, Start, End, Repeat_type)%>%
  mutate(Track = "Repeats", .after= Chromosome)

repeats$Repeat_type <- ifelse(repeats$Repeat_type == "Simple_repeat", "Simple repeat", "Others")
repeats$Repeat_type <- factor(repeats$Repeat_type, levels = c("Simple repeat", "Others"))
r_lim <- topCov + 2
r_colors <- colorRampPalette(brewer.pal(12, "Paired"))(nlevels(repeats$Repeat_type))

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
    scale_color_manual(name = "Structural variant", values = s_colors)+
    guides(color = guide_legend(order=2))+
    new_scale_color()+
  geom_point(data = loci, aes(x=Start, y = l_lim, color = Loci))+  
    scale_color_manual(name = "Loci", values = l_colors)+
    guides(color = guide_legend(order=3))+
  facet_wrap(~factor(Chromosome, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)),strip.position = "right", ncol = 2)+
  labs(y = "Normalized coverage", title = paste("Lineage:",lineage, "Sample:", sample,  sep = " "))+
  scale_y_continuous(breaks = c(1, 2)) +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 2))

ggsave("plot_simple_repeats.svg", c, height = 7, width = 13)

