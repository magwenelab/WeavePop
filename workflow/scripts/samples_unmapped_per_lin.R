log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(RColorBrewer))

#setwd("/hpc/group/magwenelab/czirion/projects/experimental_evo2/")
#metadata<- read.csv("config/sample_metadata.csv", header = TRUE)
metadata <- read.csv(snakemake@input[[1]], header = TRUE)
rownames(metadata) <- metadata$sample
metadata$group <- as.factor(metadata$group)

for (lin in levels(metadata$group)){
  #genes<-read_delim(paste(paste("results/references", lin, lin, sep ="/"), ".gff.tsv", sep = ""), col_names = TRUE, na = "N/A", show_col_types = FALSE )
  
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
    
    #file <- paste("results/samples/liftoff", samp, "unmapped_features.txt", sep = "/")
    
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
      sink(paste(snakemake@params[[4]], paste(lin, "unmapped_count.tsv", sep = "_"), sep = "/"))
      cat("There are no unmapped features in your set of samples.")
      sink()
    file.create(paste(snakemake@params[[3]], paste(lin, "unmapped.png", sep = "_"), sep = "/"))

  } else {
    
    unmapped_count <- unmapped %>%
      select(samples$sample)
    unmapped_count <- colSums(unmapped_count == 0)
    unmapped_count <-as.data.frame(unmapped_count)
    unmapped_count$sample <- rownames(unmapped_count)
    
    write.table(unmapped_count, file = paste(snakemake@params[[4]], paste(lin, "unmapped_count.tsv", sep = "_"), sep = "/"), col.names = FALSE, row.names = FALSE, quote = FALSE)

    mat <- unmapped %>%
      select(samples$sample)%>%
      mutate_all(as.integer)%>%
      as.matrix()
    
    colors <-  c( "0" = "gray", "1" = "hotpink4")
    featureCols =colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(unmapped$Feature_type)))
    names(featureCols) = unique(unmapped$Feature_type)
    split <- select(unmapped, Chromosome)
    row_ha <- rowAnnotation(Feature_type = unmapped$Feature_type, col = list(Feature_type = featureCols))

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
            show_heatmap_legend =TRUE,
            heatmap_legend_param = list(
              title = "Mapped features", at = c(0,1), 
              labels = c("Unmapped", "Mapped")))
  
    #png(paste("unmapped_", lin, ".png", sep = ""))
    svg(paste(snakemake@params[[3]], paste(lin, "unmapped.svg", sep = "_"), sep = "/"))
    draw(plot)
    dev.off()
  }
}


