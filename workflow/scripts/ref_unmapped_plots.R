log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggnewscale))

print("Reading input parameters...")
path<-snakemake@input$tsv
chromosomes_path<-snakemake@input$chroms

main_ref <- snakemake@params$main_ref
# find_repeats <- snakemake@params$find_repeats

output_path <- snakemake@output$plot
lin <- snakemake@wildcards$lineage 


print("Reading files...")
chrom_names <- read.delim(chromosomes_path, sep = ",", header = TRUE, stringsAsFactors = TRUE)
unmapped_genes <- read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

if (!is.null(snakemake@input$repeats)) {
  print("Repeats file provided...")
  print("Reading repeats file...")


} else {
  print("No repeats file provided, creating empty repeats table...")
  col_names <- c("accession", "start", "end", "repeat_type")
  repeats_table <- data.frame(matrix(ncol = length(col_names), nrow = length(chrom_names$accession)))
  colnames(repeats_table) <- col_names
  repeats_table$accession <- chrom_names$accession
  repeats_table$start <- 0
  repeats_table$end <- 1
  repeats_table$repeat_type <- "."}

print("Ordering chromosome names...")
chrom_names <- chrom_names%>%
    filter(lineage == main_ref)
unique_levels <- unique(chrom_names$chromosome)
new_order <- c(rbind(matrix(unique_levels, nrow = 2, byrow = TRUE)))
# If unique levels is an odd number, the last element is removed
if (length(unique_levels) %% 2 != 0){
  new_order <- new_order[-length(new_order)]
}
chrom_names$chromosome <- factor(chrom_names$chromosome, levels = new_order)

chroms1 <- chrom_names %>%
    filter(lineage == main_ref)%>%
    mutate(feature = paste(main_ref, "repeats", sep = " "))
chroms2 <- chrom_names %>%
    filter(lineage == main_ref)%>%
    mutate(feature = "Unmapped genes")
chroms <- rbind(chroms1, chroms2)

print("Wrangling table of genes...")
genes <- left_join(unmapped_genes, chrom_names, by = "accession")%>%
        select(chromosome,length,  start, end, gene_id, description,gene_name, feature = all_of(lin))%>%
        filter(feature == "unmapped")%>%
        mutate(feature = "Unmapped genes")

print("Counting unmapped genes...")
counts <- genes %>%
    group_by(chromosome, length)%>%
    summarize(n = n())%>%
    ungroup()%>%
    mutate(feature = "Unmapped genes")

print("Wrangling table of repeats...")
repeats<- left_join(chroms1, repeats_table, by = "accession")%>%
  select(accession, chromosome, start, end, repeat_type)

repeats$repeat_type <- ifelse(repeats$repeat_type == "Simple_repeat", "Simple repeat", "Others")
repeats$repeat_type <- factor(repeats$repeat_type, levels = c("Simple repeat", "Others"))
r_colors <- colorRampPalette(brewer.pal(12, "Paired"))(nlevels(repeats$repeat_type))
repeats$feature <- paste(main_ref, "repeats", sep = " ")

print("Plotting...")
right_limit = max(chroms$length)+50000

u <- ggplot()+
    geom_segment(data = chroms, aes(x=1, xend=length, y =feature, yend =feature),
                linewidth = 5, color = "gray95")+
    geom_segment(data = repeats, aes(x=start, xend=end, y =feature, yend =feature, color = repeat_type),
                linewidth = 5)+
    scale_color_manual(name = "Type of Repetitive\nSequences", values = r_colors)+
    geom_segment(data = genes, aes(x=start, xend=end, y =feature, yend =feature),#, color = unmapped),
                color = "black",
                linewidth = 5)+
    geom_text(data = counts, aes(x = length + 5000, y = feature, label = n, hjust = 0), size = 3)+
    scale_x_continuous(name = "Position (bp) ", labels = comma, limits = c(0, right_limit), expand = c(0.02,0))+
    facet_wrap(~chromosome, ncol=2, strip.position = "right")+
    theme_classic()+
    theme(legend.position = "bottom")+
    labs(title = paste("Location of Genes of the Main Reference Unmapped in ", lin),
        y = "")

print("Saving plot...")
ggsave(output_path, u, width = 10, height = 6, dpi = 600)

print("Done!")