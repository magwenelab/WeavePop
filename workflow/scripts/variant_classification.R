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
lineage_barplot_path = snakemake@output$plot
# variant_classification_path = "/FastData/czirion/WeavePop/test_bigger/results/04.Intermediate_files/02.Dataset/snpeff/VNI_variants_classification.tsv"
# lineage_barplot_path = paste("./", lineage, "_variants_summary_barplot.png", sep = "")

print("Reading files...")
variants = read.delim(variants_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
effects = read.delim(effects_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
presence = read.delim(presence_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
metadata = read.delim(metadata_path, sep = ",", header = TRUE, stringsAsFactors = TRUE)
chromosomes = read.delim(chromosomes_path, sep = ",", header = TRUE, stringsAsFactors = TRUE)

# Assign one impact per variant in order of priority HIGH > MODERATE > LOW > MODIFIEReffs <- effects %>%
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

# Add the impact table to the variants table and add status (private to reference or sample or Non-private in multiple samples) to variants table
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

vars_type <- vars_classification%>%
    left_join(chromosomes, by = c("accession", "lineage"))%>% 
    select(var_id, lineage, chromosome, impact, status )
total_vars = nrow(vars_type)

# Barplot summary of number of variants per status and impact in lineage VERSION 1
b <- ggplot(vars_type, aes(x = status, fill = impact))+
    geom_bar(stat = "count", position = "dodge")+
    geom_text(stat= "count", position = position_dodge(width =0.9), 
                aes(label = scales::comma(after_stat(count))),
                vjust = -0.5,
                size = 3)+
    scale_y_continuous(name = "Number of variants", labels = comma)+
    labs(title = paste("Number of variants in each status and impact of lineage ", lineage),
         subtitle = paste("Total number of variants:", scales::comma(total_vars)),
         x = "", 
         fill = "Impact")+
    theme_classic()

# Boxplot summary of number of variants per sample per status and impact VERSION 2
vars_type_sample <- left_join(presence, vars_type, by = "var_id")%>%
    group_by(sample, lineage, status, impact)%>%
    summarize(n_variants = n())

vars_type_impact <-  vars_type %>%
    group_by(status, impact)%>%
    summarize(n_variants = n())

g <- ggplot(vars_type_sample, aes(x = impact, y = n_variants))+
    expand_limits(y = 0)+
    geom_quasirandom(aes(color = impact), alpha = 0.5)+
    geom_boxplot(alpha = 0)+
    facet_wrap(~status, nrow = 1, scales = "free")+
    scale_y_continuous(name = "Number of variants per sample", labels = comma)+
    labs(title = paste("Number of variants per sample in each status and impact of lineage ", lineage),
         subtitle = paste("Total number of variants:", scales::comma(total_vars)),
         x = "", 
         color = "Impact")+
    theme_classic()

# Boxplot and barplot summaries of number of variants per sample per status and impact VERSION 3
## It would be better to just add the number of the barplot to the top of the boxplots
c <- g / b

print("Saving variant classification table...")
write_delim(vars_classification, variant_classification_path, delim = "\t")

print("Saving plot...")
ggsave(lineage_barplot_path, c, width = 8, height = 5) # Choose g, b, or c
print("Done!")