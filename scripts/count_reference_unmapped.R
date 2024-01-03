log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))
library(RColorBrewer)

# genes<-read_delim("references/FungiDB-65_CneoformansH99.gff.tsv", col_names = TRUE, na = "N/A", show_col_types = FALSE )
genes<-read_delim(snakemake@input[[1]], col_names = TRUE, na = "N/A", show_col_types = FALSE )
genes<- genes %>% 
  # filter(str_detect(primary_tag, "gene" ))%>%
  as.data.frame()
rownames(genes)<- genes$ID
# lins <- read.csv("files/lineage_references.csv", header = TRUE)
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
  select(Chromosome = seq_id, Feature_type = primary_tag, lins$group)%>%
  filter(rowSums(. == 0) > 0)

unmapped_count <- unmapped %>%
  select(lins$group)
unmapped_count <- colSums(unmapped_count == 0)
unmapped_count <-as.data.frame(unmapped_count)
unmapped_count$lineage <- rownames(unmapped_count)

# write.table(unmapped_count, file = "unmapped_count.tsv", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(unmapped_count, file = snakemake@output[[1]],  col.names = FALSE, row.names = FALSE, quote = FALSE)

mat <- unmapped %>%
  select(lins$group)%>%
  mutate_all(as.integer)%>%
  as.matrix()
  
colors <-  c( "0" = "gray", "1" = "black")
featureCols =colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(unmapped$Feature_type)))
names(featureCols) = unique(unmapped$Feature_type)
split <- select(unmapped, Chromosome)
row_ha <- rowAnnotation(Feature_type = unmapped$Feature_type, col = list(Feature_type = featureCols))

# svg("unmapped.svg",width=16,height=25)
svg(snakemake@output[[2]],width=16,height=25)
Heatmap(mat, 
        name = "Mapped features",
        col = colors,
        show_row_names = TRUE,
        cluster_rows = TRUE,
        row_split = split, 
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        row_title_rot = 0,
        right_annotation = row_ha,
        row_names_gp = gpar(fontsize = 5))
dev.off()


