log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(patchwork))

print("Reading input parameters...")
samp = snakemake@wildcards$sample
window_size = snakemake@params$window_size
variants_path = snakemake@input$variants
classif_path = snakemake@input$classif
presence_path = snakemake@input$presence
chromosomes_path = snakemake@input$chromosomes
metadata_input <- snakemake@input$metadata
summary_path = snakemake@output$summary
category_path = snakemake@output$category
impact_path = snakemake@output$impact


print("Reading files...")
variants = read.delim(variants_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
classif = read.delim(classif_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
presence = read.delim(presence_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
chromosomes = read.delim(chromosomes_path, sep = ",", header = TRUE)
metadata <- read.delim(metadata_input, sep = ",", header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))

classif$category <- factor(classif$category, levels = c("Reference private","Private", "Non-private"))
chromosomes$chromosome  <- factor(chromosomes$chromosome, levels = unique(chromosomes$chromosome))

print("Obtaining lineage of sample...")


lineage_name <- as.character(metadata$lineage[metadata$sample == samp])
strain_name <- as.character(metadata$strain[metadata$sample == samp])

print("Merge variant description and classification")
variants <- left_join(variants, classif, by = "var_id")

print("Add chromosome name and creating windows")
vars_windows <- variants %>%
    left_join(chromosomes, by = c("accession", "lineage"))%>%
    group_by(chromosome) %>%
    mutate(window = floor(pos / window_size) + 1)%>%
    ungroup()%>%
    group_by(chromosome, window) %>%
    mutate(window_start = (window - 1) * window_size + 1,
            window_end = window * window_size)%>%
    ungroup()
    
print("Filtering variants of sample")
vars_sample <- presence %>%
    filter(sample == samp)%>%
    droplevels()

n_vars <- length(unique(vars_sample$var_id))
vars <- vars_windows %>%
    filter(var_id %in% vars_sample$var_id)

if (length(unique(presence$sample)) > 1){
    print("Removing variants private to reference genome")
    vars <- vars%>%
        filter(category != "Reference private")%>%
        droplevels()
} else {
    print("Only one sample in lineage, plotting all variants")
}

print("Getting number of variants per window for each category")
category_density <- vars %>% 
    group_by(chromosome, window,window_start, window_end, category, length)%>%
    summarize(n = n())

print("Getting number of variants per window for each impact")
effs_density <- vars%>% 
    group_by(chromosome, window,window_start, window_end, impact, length)%>%
    summarize(n = n())

print("Plotting category density")
p <- ggplot(category_density)+
    geom_segment(aes(x=1, xend=length, y=-1, yend=-1), linewidth = 0.1, color = "black")+
    geom_col(aes(x = window_start, y = n, color = category, fill = category))+
    facet_grid(chromosome~category)+
    scale_x_continuous(name = "Position (bp) ", labels = comma)+
    labs(title = "Number of Variants in Windows Along Chromosomes by Privateness Category",
        subtitle= paste("Sample:", samp, "Strain:", strain_name, " Lineage:", lineage_name, " Window size:", window_size, " Total variants: ", scales::comma(n_vars)),
            y = "Number of variants")+
    theme_classic()+
    theme(legend.position = "none")

print("Plotting impact density")
e <- ggplot(effs_density)+
    geom_segment(aes(x=1, xend=length, y=-1, yend=-1), linewidth = 0.1, color = "black")+
    geom_col(aes(x = window_start, y = n, color = impact, fill = impact))+
    facet_grid(chromosome~impact, scales = "free_x")+#, ncol =1, strip.position = "right")+
    scale_x_continuous(name = "Position (bp) ", labels = comma)+
    labs(title="Number of Variants in Windows Along Chromosomes by Impact",
        subtitle= paste("Sample:", samp, "Strain:", strain_name, " Lineage:", lineage_name, " Window size:", window_size," Total variants: ", scales::comma(n_vars)),
            y = "Number of variants")+
    theme_classic()+
    theme(legend.position = "none")

print("Plotting summary of number of variants")
total_vars = nrow(vars)
b <- ggplot(vars, aes(x = category, fill = impact))+
    geom_bar(stat = "count", position = "dodge")+
    geom_text(stat= "count", position = position_dodge(width =0.9), 
                aes(label = scales::comma(after_stat(count))),
                vjust = -0.5,
                size = 3)+
    scale_y_continuous(name = "Number of variants", labels = comma)+
    labs(title = "Number of Variants per Impact and Privateness Category",
        subtitle= paste("Sample:", samp, "Strain:", strain_name, " Lineage:", lineage_name, " Total variants:", scales::comma(total_vars)),
        x = "", 
        fill = "Impact")+
    theme_classic()

print("Saving plots...")
ggsave(summary_path, b, width = 8, height = 5)
ggsave(category_path, p, width = 16, height = 9)
ggsave(impact_path, e, width = 16, height = 9)
print("Done!")