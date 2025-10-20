log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(patchwork))

print("Reading input parameters...")
effects_path  = snakemake@input$effects
presence_path = snakemake@input$presence
variants_path = snakemake@input$variants
metadata_path = snakemake@input$metadata
chromosomes_path = snakemake@input$chromosomes
lineage = snakemake@wildcards$lineage

# presence_path = "/FastData/czirion/WeavePop/test_bigger/results/04.Intermediate_files/02.Dataset/snpeff/VNI_presence.tsv"
# variants_path = "/FastData/czirion/WeavePop/test_bigger/results/04.Intermediate_files/02.Dataset/snpeff/VNI_variants.tsv"
# metadata_path = "/FastData/czirion/WeavePop/test_bigger/results/02.Dataset/metadata.csv"
# chromosomes_path = "/FastData/czirion/WeavePop/test_bigger/results/02.Dataset/chromosomes.csv"
# lineage = "VNBI"

variant_classification_path = snakemake@output$tsv

# variant_classification_path = "/FastData/czirion/WeavePop/test_bigger/results/04.Intermediate_files/02.Dataset/snpeff/VNI_variants_classification.tsv"
# lineage_barplot_path = paste("./", lineage, "_variants_summary_barplot.png", sep = "")

print("Reading files...")
variants = read.delim(variants_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
effects = read.delim(effects_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
presence = read.delim(presence_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
metadata = read.delim(metadata_path, sep = ",", header = TRUE, stringsAsFactors = TRUE)
chromosomes = read.delim(chromosomes_path, sep = ",", header = TRUE, stringsAsFactors = TRUE)

print("Assigning one impact per variant in order of priority HIGH > MODERATE > LOW > MODIFIER")
effs <- effects %>%
    group_by(var_id, impact)%>%
    summarize(n_effects = n())%>%
    mutate(one_effect = 1)%>%
    select(-n_effects)%>%
    pivot_wider(names_from = impact, values_from = one_effect, values_fill = 0) %>%
    mutate(impact = ifelse(HIGH == 1, "High",
                         ifelse(MODERATE == 1, "Moderate",
                             ifelse(LOW == 1, "Low", "Modifier"))))%>%
    select(var_id, impact)%>%
    ungroup()
effs$impact <- factor(effs$impact, levels = c("High", "Moderate", "Low", "Modifier"))    

print("Adding impact and status (private to reference or sample or Non-private in multiple samples)")
samples_in_lineage <- metadata %>%
    group_by(lineage)%>%
    summarize(samples_in_lineage = length(unique(sample)))

presence_summary <- presence %>%
    group_by(var_id)%>%
    summarize(n_samples = n())

vars_classification <- left_join(variants, effs, by = "var_id") %>%
    left_join(samples_in_lineage, by = "lineage") %>%
    left_join(presence_summary, by = "var_id") %>%
    mutate(status = ifelse(n_samples == samples_in_lineage, "Reference private", ifelse(n_samples == 1, "Private", "Non-private") ))%>%
    select(-c(samples_in_lineage, n_samples))
vars_classification$status <- factor(vars_classification$status, levels = c("Reference private","Private", "Non-private"))

print("Saving variant classification table...")
write_delim(vars_classification, variant_classification_path, delim = "\t")
