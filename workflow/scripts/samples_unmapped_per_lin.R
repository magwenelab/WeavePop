log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(ComplexHeatmap)
library(RColorBrewer)
#setwd("/hpc/group/magwenelab/czirion/projects/experimental_evo/")
#metadata<- read.csv("files/sample_metadata.csv", header = TRUE)
metadata <- read.csv(snakemake@input[[1]], header = TRUE)
rownames(metadata) <- metadata$sample
metadata$group <- as.factor(metadata$group)

for (lin in levels(metadata$group)){
    REFDIR / "{lineage}" / "{lineage}.gff"
  genes<-read_delim(paste(paste(snakemake@params[[1]], lin, lin, sep ="/"), ".gff.tsv", sep = ""), col_names = TRUE, na = "N/A", show_col_types = FALSE )
  
  genes<- genes %>% 
    select(Chromosome = seq_id, Feature_type = primary_tag, ID, description = matches("description|product"), contains("Name"))%>%
    filter(str_detect(Feature_type, "gene" ))%>%
    as.data.frame()
  toFeature <- colnames(genes)[! colnames(genes) %in% c("Chromosome", "Feature_type")]
  genes <- unite(genes, Feature, all_of(toFeature), sep = " ", remove = FALSE, na.rm = TRUE)
  rownames(genes)<- genes$Feature 

  samples <- metadata %>%
    filter(group == lin)

  for (samp in samples$sample){
    file <- paste(snakemake@params[[2]], samp, "unmapped_features.txt", sep = "/")
    df<- read.csv(file, header = FALSE, col.names = c("ID"), colClasses = "character")
    genes <- genes %>%
      mutate(!!samp := ifelse(ID %in% df$ID, 0, 1))
  }
  
  unmapped <- genes %>% 
    select(Chromosome, Feature_type, samples$sample)%>%
    filter(rowSums(. == 0) > 0)
  
  if(nrow(unmapped)== 0){
    print('There are no unmapped features in your set of samples.')
    file.create(snakemake@output[[1]])
    file.create(snakemake@output[[2]])
  } else {
    
    unmapped_count <- unmapped %>%
      select(samples$sample)
    unmapped_count <- colSums(unmapped_count == 0)
    unmapped_count <-as.data.frame(unmapped_count)
    unmapped_count$sample <- rownames(unmapped_count)
    
    write.table(unmapped_count, file = "unmapped_count.tsv", col.names = FALSE, row.names = FALSE, quote = FALSE)
    #write.table(unmapped_count, file = snakemake@output[[1]],  col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    mat <- unmapped %>%
      select(samples$sample)%>%
      mutate_all(as.integer)%>%
      as.matrix()
    
    colors <-  c( "0" = "gray", "1" = "black")
    featureCols =colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(unmapped$Feature_type)))
    names(featureCols) = unique(unmapped$Feature_type)
    split <- select(unmapped, Chromosome)
    row_ha <- rowAnnotation(Feature_type = unmapped$Feature_type, col = list(Feature_type = featureCols))
    pwidth = 5 + 0.5 * nlevels(as.factor(metadata$sample))
    pheight = 3 + 0.05 * nrow(unmapped)

    plot <- Heatmap(mat, 
            column_title = lin,
            col = colors,
            show_row_names = TRUE,
            cluster_rows = TRUE,
            row_split = split, 
            show_row_dend = FALSE,
            show_column_dend = FALSE,
            row_title_rot = 0,
            right_annotation = row_ha,
            row_names_gp = gpar(fontsize = 5),
            column_names_gp = gpar(fontsize = 5),
            show_heatmap_legend = FALSE)

   png(paste(snakemake@params[[1]], "unmapped_", lin, ".png", sep = ""))
   draw(plot)
   dev.off()
  }
}

