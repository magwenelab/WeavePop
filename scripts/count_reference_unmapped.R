log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
library(pheatmap)
library(RColorBrewer)

#genes<-read_delim("references/South.gff.tsv", col_names = TRUE, na = "N/A", show_col_types = FALSE )
genes<-read_delim(snakemake@input[[1]], col_names = TRUE, na = "N/A", show_col_types = FALSE )
genes<- genes %>% 
  filter(str_detect(primary_tag, "gene" ))%>%
  as.data.frame()
rownames(genes)<- genes$ID
#lins <- read.csv("files/lineage_references.csv", header = TRUE)
lins <- read.csv(snakemake@input[[2]], header = TRUE)

# for (lin in lins$group){
#   file <- paste("references/", lin, "_unmapped_features.txt", sep = "")
#   df<- read.csv(file, header = FALSE, col.names = c("ID"), colClasses = "character")
#   genes <- genes %>%
#     mutate(!!lin := ifelse(ID %in% df$ID, 0, 1))
# }
for (lin in lins$group){
  file <- paste(snakemake@config[["reference_directory"]], lin, "_unmapped_features.txt", sep = "")
  df<- read.csv(file, header = FALSE, col.names = c("ID"), colClasses = "character")
  genes <- genes %>%
    mutate(!!lin := ifelse(ID %in% df$ID, 0, 1))
}

unmapped <- genes %>% 
  select(lins$group)%>%
  filter(rowSums(. == 0) > 0)

unmapped_count <- colSums(unmapped == 0)
write.table(unmapped_count, file = snakemake@output[[1]], col.names = FALSE, row.names = TRUE)

unmapped_matrix<-as.matrix(unmapped)

genes <- genes %>% 
  filter(ID %in% rownames(unmapped))%>%
  select(Feature_type = primary_tag, Chromosome = seq_id)

lg.brks<-c(0,1) 
lb.brks<-c("Unmapped","Mapped") 
plot <- pheatmap(unmapped_matrix,
                 color =  c("gray70", "gray20"),
                 annotation_row = genes,
                 legend_breaks = lg.brks, 
                 legend_labels = lb.brks,
                 cluster_cols = FALSE,
                 cluster_rows = FALSE,
                 cellwidth = 10,
                 cellheight = 2,
                 fontsize_row= 2,
                 treeheight_row = 0,
                 treeheight_col = 0)
ggsave(snakemake@output[[2]], plot, dpi = 50, units = "cm", height = 23, width = 16)

# Missing to add gene name labels