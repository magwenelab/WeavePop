log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(pheatmap)

#genes<-read_delim("references/FungiDB-65_CneoformansH99.gff.tsv", col_names = TRUE, na = "N/A", show_col_types = FALSE )
genes<-read_delim(snakemake@input[[1]], col_names = TRUE, na = "N/A", show_col_types = FALSE )
genes<- genes %>% 
  filter(str_detect(primary_tag, "gene" ))%>%
  as.data.frame()
rownames(genes)<- genes$ID

#samples <- read.csv("files/sample_metadata.csv", header = TRUE)
samples <- read.csv(snakemake@input[[2]], header = TRUE)
rownames(samples) <- samples$sample
for (samp in samples$sample){
  file <- paste("genomes-annotations/", samp, "/unmapped_features.txt", sep = "")
  df<- read.csv(file, header = FALSE, col.names = c("ID"), colClasses = "character")
  genes <- genes %>%
    mutate(!!samp := ifelse(ID %in% df$ID, 0, 1))
}

unmapped <- genes %>% 
  select(samples$sample)%>%
  filter(rowSums(. == 0) > 0)

unmapped_count <- colSums(unmapped == 0)
#write.table(unmapped_count, file = "results/unmapped_count.txt", col.names = FALSE, row.names = TRUE)
write.table(unmapped_count, file = snakemake@output[[1]], col.names = FALSE, row.names = TRUE)

unmapped_matrix<-as.matrix(unmapped)
cols_to_remove <- which(colSums(unmapped_matrix == 1) == nrow(unmapped_matrix))
unmapped_matrix <- unmapped_matrix[, -cols_to_remove]

samples <- samples %>% 
  filter(sample %in% colnames(unmapped_matrix))%>%
  select(Lineage = group)

genes <- genes %>% 
  filter(ID %in% rownames(unmapped))%>%
  select(Feature_type = primary_tag, Chromosome = seq_id)

lg.brks<-c(0,1) 
lb.brks<-c("Unmapped","Mapped") 
plot <- pheatmap(unmapped_matrix,
                 color =  c("gray70", "gray20"),
                 annotation_row = genes,
                 annotation_col = samples,
                 legend_breaks = lg.brks, 
                 legend_labels = lb.brks,
                 cluster_cols = TRUE,
                 cluster_rows = FALSE,
                 cellwidth = 10,
                 cellheight = 2,
                 fontsize_row= 2,
                 treeheight_row = 0,
                 treeheight_col = 0)

ggsave(snakemake@output[[2]], plot, units = "cm", height = 20, width = 50)
#ggsave("results/unmapped.svg", plot, units = "cm", height = 20, width = 50)
