log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(svglite))

# samples <- read.csv("config/sample_metadata.csv", header = TRUE)
samples <- read.csv(snakemake@input[[1]], header = TRUE)
rownames(samples) <- samples$sample
# genes<-read_delim("results/references/FungiDB-65_CneoformansH99.tsv", col_names = TRUE, na = "N/A", show_col_types = FALSE )
genes<-read_delim(snakemake@input[[2]], col_names = TRUE, na = "N/A", show_col_types = FALSE )
genes<- genes %>% 
  select(Chromosome = seq_id, Feature_type = primary_tag, ID, description = matches("description|product"), contains("Name"))%>%
  filter(str_detect(Feature_type, "gene" ))%>%
  as.data.frame()
toFeature <- colnames(genes)[! colnames(genes) %in% c("Chromosome", "Feature_type")]
genes <- unite(genes, Feature, all_of(toFeature), sep = " ", remove = FALSE, na.rm = TRUE)
rownames(genes)<- genes$Feature
# for (samp in samples$sample){
#   file <- paste("results/samples/liftoff", samp, "unmapped_features.txt", sep = "/")
#   df<- read.csv(file, header = FALSE, col.names = c("ID"), colClasses = "character")
#   genes <- genes %>%
#     mutate(!!samp := ifelse(ID %in% df$ID, 0, 1))
# }

for (samp in samples$sample){
  file <- paste(snakemake@params[[1]], samp, "unmapped_features.txt", sep = "/")
  df<- read.csv(file, header = FALSE, col.names = c("ID"), colClasses = "character")
  genes <- genes %>%
    mutate(!!samp := ifelse(ID %in% df$ID, 0, 1))
}
unmapped <- genes %>% 
  select(Chromosome, Feature_type, samples$sample)%>%
  filter(rowSums(. == 0) > 0)

if(nrow(unmapped)== 0){
  print('There are no unmapped features in your set of samples.')
  sink(snakemake@output[[1]])
  cat("There are no unmapped features in your set of samples.")
  sink()
  file.create(snakemake@output[[2]])
} else {
  unmapped_count <- unmapped %>%
    select(samples$sample)
  unmapped_count <- colSums(unmapped_count == 0)
  unmapped_count <-as.data.frame(unmapped_count)
  unmapped_count$sample <- rownames(unmapped_count)

  #write.table(unmapped_count, file = "unmapped_count.tsv", col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(unmapped_count, file = snakemake@output[[1]],  col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  mat <- unmapped %>%
    select(samples$sample)%>%
    mutate_all(as.integer)%>%
    as.matrix()

  colors <-  c( "0" = "gray", "1" = "hotpink4")
  featureCols =colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(unmapped$Feature_type)))
  names(featureCols) = unique(unmapped$Feature_type)
  linCols =colorRampPalette(brewer.pal(12, "Paired"))(length(unique(samples$lineage)))
  names(linCols) = unique(samples$lineage)
  row_split <- select(unmapped, Chromosome)
  col_split <- select(samples, lineage)
  row_ha <- rowAnnotation(Feature_type = unmapped$Feature_type, col = list(Feature_type = featureCols))
  #col_ha <- HeatmapAnnotation(Lineage = samples$lineage , col = list(Lineage = linCols))

  plot <-   Heatmap(mat, 
          col = colors,
          show_row_names = TRUE,
          cluster_rows = TRUE,
          row_split = row_split, 
          column_split = col_split,
          show_row_dend = FALSE,
          show_column_dend = FALSE,
          row_title_rot = 0,
          column_title_rot = 90,
          right_annotation = row_ha,
          #top_annotation = col_ha,
          row_names_gp = gpar(fontsize = 5),
          column_names_gp = gpar(fontsize = 3),
          show_heatmap_legend = TRUE,
          heatmap_legend_param = list(
              title = "Mapped features", at = c(0,1), 
              labels = c("Unmapped", "Mapped")))

  #png("unmapped.png",width=pwidth,height=pheight)
  svg(snakemake@output[[2]], width = 13.3, height= 7.5)
  draw(plot)
  dev.off()

}
