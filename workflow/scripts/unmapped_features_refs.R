log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))
library(RColorBrewer)

# lins <- read.csv("config/sample_metadata.csv", header = TRUE)
lins <- read.csv(snakemake@input[[1]], header = TRUE)

# genes<-read_delim("results/references/FungiDB-65_CneoformansH99.tsv", col_names = TRUE, na = c("NA","N/A", ""), show_col_types = FALSE )
genes<-read_delim(snakemake@input[[2]], col_names = TRUE, na = c("NA","N/A", ""), show_col_types = FALSE )
genes<- genes %>% 
  select(Chromosome = seq_id, Feature_type = primary_tag, ID, description = matches("description|product"), contains("Name"))%>%
  filter(str_detect(Feature_type, "gene" ))%>%
  as.data.frame()
toFeature <- colnames(genes)[! colnames(genes) %in% c("Chromosome", "Feature_type")]
genes <- unite(genes, Feature, all_of(toFeature), sep = " ", remove = FALSE, na.rm = TRUE)
rownames(genes)<- genes$Feature


# for (lin in levels(as.factor(lins$group))){
#   file <- paste("results/references/", lin, "/unmapped_features.txt", sep = "")
#   df<- read.csv(file, header = FALSE, col.names = c("ID"), colClasses = "character")
#   genes <- genes %>%
#     mutate(!!lin := ifelse(ID %in% df$ID, 0, 1))
# }
for (lin in levels(as.factor(lins$group))){
 file <- paste(snakemake@params[[1]], lin, "unmapped_features.txt", sep = "/")
 df<- read.csv(file, header = FALSE, col.names = c("ID"), colClasses = "character")
 genes <- genes %>%
   mutate(!!lin := ifelse(ID %in% df$ID, 0, 1))
}

unmapped <- genes %>% 
  select(Chromosome, Feature_type, lins$group)%>%
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


pwidth = 5 + 0.5 * nlevels(as.factor(lins$group))
pheight = 3 + 0.05 * nrow(unmapped)
#png("unmapped.png",width=pwidth,height=pheight)
png(snakemake@output[[2]],width=pwidth,height=pheight)
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
        row_names_gp = gpar(fontsize = 5),
        show_heatmap_legend = FALSE)
dev.off()