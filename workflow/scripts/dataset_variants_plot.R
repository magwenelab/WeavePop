log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(RColorBrewer))


print("Reading input parameters...")

metadata_path <- snakemake@input$metadata
variant_classification_path <- snakemake@input$classif
variants_path <- snakemake@input$variants
presence_path <- snakemake@input$presence

output <- snakemake@output$plot

# variant_classification_path <- "/FastData/czirion/WeavePop/test/results/02.Dataset/snpeff/variant_classification.tsv"
# variants_path <- "/FastData/czirion/WeavePop/test/results/02.Dataset/snpeff/variants.tsv"
# presence_path <- "/FastData/czirion/WeavePop/test/results/02.Dataset/snpeff/presence.tsv"

# output_path <- "plot.png"

print("Reading files...")
metadata <- read.delim(metadata_path, sep = ",", header = TRUE, stringsAsFactors = TRUE)
vars_classification <- read.delim(variant_classification_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
variants <- read.delim(variants_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
presence <- read.delim(presence_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)

if ("strain" %in% colnames(metadata)){
    metadata <- metadata %>%
        select(sample, strain, ref_genome) %>%
        mutate(name = paste(strain, sample, sep = " ")) 
    } else {
    metadata <- metadata %>%
        select(sample, ref_genome) %>%
        mutate(name = sample)}
        
metadata <- metadata %>%
    arrange(name)%>%
    mutate(name = factor(name, levels = unique(name)))

vars_classification$category <- factor(vars_classification$category, levels = c("Reference private","Private", "Non-private"))
vars_classification$impact<- factor(vars_classification$impact, levels = c("High","Moderate","Low", "Modifier"))

vars_classification <- left_join(variants, vars_classification, by = "var_id")%>%
    filter(category != "Reference private")



paired <- brewer.pal(12, "Paired")
i_colors <- c(paired[c(6,8,2,10)])
names(i_colors) <- levels(vars_classification$impact)
c_colors <- c(paired[c(4,1)])
names(c_colors) <- c("Private",  "Non-private")

vars_type <- vars_classification%>%
    select(var_id, ref_genome, impact, category )

vars_type_sample <- left_join(vars_type, presence, by = "var_id")%>%
    left_join(metadata, by = c("sample", "ref_genome"))%>%
    group_by(sample, name, ref_genome, category, impact)%>%
    summarize(n_variants = n())%>%
    ungroup()

vars_category_sample <- vars_type_sample %>%
    group_by(sample,name, category, ref_genome)%>%
    summarize(total_category = sum(n_variants))%>%
    ungroup()%>%
    group_by(sample)%>%
    mutate(percent_category = total_category / sum(total_category) * 100)

vars_impact_sample <- vars_type_sample %>%
    group_by(sample,name, impact, ref_genome)%>%
    summarize(total_impact = sum(n_variants))%>%
    ungroup()%>%
    group_by(sample)%>%
    mutate(percent_impact = total_impact / sum(total_impact) * 100)

vars_total_sample <- vars_type_sample %>%
    group_by(sample,name, ref_genome)%>%
    summarize(n_variants = sum(n_variants))%>%
    ungroup()

t <-ggplot() +
        geom_point(data = vars_total_sample, aes(x = name, y = n_variants)) +
        facet_grid(~ ref_genome, scales = "free", space = "free_x") +
        theme_bw()+
        theme(panel.background = element_blank(), 
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())+
        labs(x = "", y = "# Variants",  title = "Number of Filtered Variants")+
        scale_y_continuous(labels = scales::comma)

c <-ggplot() +
        geom_bar(data = vars_category_sample, aes(x = name, y = percent_category, fill = category), stat = "identity") +
        facet_grid(~ ref_genome, scales = "free", space = "free_x") +
        theme(panel.background = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            strip.text = element_blank(),
            panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())+
        labs(x = "", y = "% Variants", fill = "Metric", title = "Percentage of Variants by Category")+
        scale_fill_manual(values = c_colors, name = "")



i <-ggplot() +
        geom_bar(data = vars_impact_sample, aes(x = name, y = percent_impact, fill = impact), stat = "identity") +
        facet_grid(~ ref_genome, scales = "free", space = "free_x") +
        theme(panel.background = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            strip.text = element_blank(),
            panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 1),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5))+
        labs(x = "", y = "% Variants", fill = "Metric", title = "Percentage of Variants by Impact")+
        scale_fill_manual(values = i_colors, name = "")
plot <- t/c/i

print("Saving plot...")
n_samples <- nrow(metadata)

if (n_samples <= 200) {
    gscale <- 1
    gwidth <- 6 + n_samples * 0.03
} else if (n_samples > 200 & n_samples <= 400) {
    gwidth <- 16
    gscale <- 2
} else if ( n_samples > 400 & n_samples <= 1000) {
    gwidth <- 18
    gscale <- 2.5
} else {
    gwidth <- 20
    gscale <- 3
}

gheight <- 7

ggsave(output, plot = plot, units = "in", height = gheight, width = gwidth, scale = gscale, limitsize = FALSE)

print("Done!")