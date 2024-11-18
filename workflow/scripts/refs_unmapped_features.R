log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries...")
suppressPackageStartupMessages(library(tidyverse))

print("Reading files...")
metadata <- read.csv(snakemake@input[[1]], header = TRUE)
lineages_path <- snakemake@params[[1]]
genes<-read_delim(snakemake@input[[2]], col_names = TRUE, na = c("NA","N/A", ""), show_col_types = FALSE )

print("Getting lineages...")
lineages <- unique(metadata$lineage)

print("Selecting columns of main reference...")
genes<- genes %>% 
  select(accession= seq_id, start, end, primary_tag, gene_id = ID, description , contains("Name"))%>%
  rename(gene_name = Name)

print("Reading list of unmapped features of each lineage and joining to annotation table...")
for (lin in lineages){
  file <- paste("/FastData/czirion/DiversityPipeline/results_1/04.Intermediate_files/03.References/", lin, "/unmapped_features.txt", sep = "")
  df<- read.csv(file, header = FALSE, col.names = c("gene_id"), colClasses = "character")
  df[lin]<- "unmapped"
  genes <- genes %>%
    left_join(df, by = "gene_id")
}

print("Filtering only unmapped features...")
unmapped <- genes %>%
  filter_at(vars(all_of(lineages)), any_vars(!is.na(.)))

print("Saving output...")
write_tsv(unmapped, file = snakemake@output[[1]], col_names = TRUE)
print("Done!")
