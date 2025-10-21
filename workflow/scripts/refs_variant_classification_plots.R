log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(patchwork))

print("Reading input parameters...")

chromosomes_path = snakemake@input$chromosomes #"/FastData/czirion/WeavePop_Cneoformans/Crypto_Desjardins/results/02.Dataset/chromosomes.csv"
variant_classification_path = snakemake@input$variants#"/FastData/czirion/WeavePop_Cneoformans/Crypto_Desjardins/results/04.Intermediate_files/02.Dataset/snpeff/VNBI_variant_classification.tsv"#
presence_path = snakemake@input$presence#"/FastData/czirion/WeavePop_Cneoformans/Crypto_Desjardins/results/04.Intermediate_files/02.Dataset/snpeff/VNBI_presence.tsv"#

lineage=snakemake@wildcards$lineage#"VNBI" #
window_size =snakemake@params$window_size#500 # 

lineage_barplot_path = snakemake@output$plot#"/FastData/czirion/WeavePop_Cneoformans/Crypto_Desjardins/results/02.Dataset/plots/VNBI_variant_summary.png"#
ref_private_path = snakemake@output$plot_density#"/FastData/czirion/WeavePop_Cneoformans/Crypto_Desjardins/results/02.Dataset/plots/VNBI_reference_variants.png"#
   

print("Reading files...")
vars_classification = read.delim(variant_classification_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
presence = read.delim(presence_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
chromosomes = read.delim(chromosomes_path, sep = ",", header = TRUE, stringsAsFactors = TRUE)

vars_classification$status <- factor(vars_classification$status, levels = c("Reference private","Private", "Non-private"))
vars_classification$impact<- factor(vars_classification$impact, levels = c("High","Moderate","Low", "Modifier"))

print("Ordering chromosome names...")

unique_levels <- unique(chromosomes$chromosome)
new_order <- c(rbind(matrix(unique_levels, nrow = 2, byrow = TRUE)))
# If unique levels is an odd number, the last element is removed
if (length(unique_levels) %% 2 != 0){
  new_order <- new_order[-length(new_order)]
}
chromosomes$chromosome <- factor(chromosomes$chromosome, levels = new_order)

print("Creating summary plot...")
vars_type <- vars_classification%>%
    left_join(chromosomes, by = c("accession", "lineage"))%>% 
    select(var_id, lineage, chromosome, impact, status )

total_vars = nrow(vars_type)

vars_type_sample <- left_join(presence, vars_type, by = "var_id")%>%
    group_by(sample, lineage, status, impact)%>%
    summarize(n_variants = n())%>%
    ungroup()

vars_type_sample_max <- vars_type_sample %>%
    group_by(status)%>%
    summarize(max_variants = max(n_variants) + max(n_variants)/10)

vars_summary <- vars_type %>%
    group_by(impact, status)%>%
    summarize(n_variants = n())%>%
    left_join(vars_type_sample_max, by = "status")

g <- ggplot(vars_type_sample, aes(x = impact, y = n_variants))+
    expand_limits(y = 0)+
    geom_quasirandom(aes(color = impact), alpha = 0.5)+
    geom_boxplot(alpha = 0)+
    geom_text(data = vars_summary, aes(y = max_variants, x = impact, label = scales::comma(n_variants)))+
    facet_wrap(~status, nrow = 1, scales = "free")+
    scale_y_continuous(name = "Number of variants per sample", labels = comma)+
    labs(title = paste("Number of variants per sample in each status and impact of lineage ", lineage),
         subtitle = paste("Total number of variants:", scales::comma(total_vars)),
         x = "", 
         color = "Impact")+
    theme_classic()

print("Making density plot of reference private variants")

vars_windows <- vars_classification %>%
    left_join(chromosomes, by = c("accession", "lineage"))%>%
    group_by(chromosome) %>%
    mutate(window = floor(pos / window_size) + 1)%>%
    ungroup()%>%
    group_by(chromosome, window) %>%
    mutate(window_start = (window - 1) * window_size + 1,
            window_end = window * window_size)%>%
    ungroup()%>%
    filter(status == "Reference private")

status_density <- vars_windows %>% 
    group_by(chromosome, window,window_start, window_end, status, length)%>%
    summarize(n = n())

n_vars <- length(unique(vars_windows$var_id))

print("Plotting status density")
p <- ggplot()+
    geom_segment(data = chromosomes, 
                aes(x=1, xend=length, y=-1, yend=-1), 
                linewidth = 0.1, color = "black")+
    geom_col(data = status_density, 
            aes(x = window_start, y = n), 
            color = "black", fill = "black")+
    facet_wrap(~chromosome, ncol = 2, strip.position = "right")+
    scale_x_continuous(name = "Position (bp) ", labels = comma)+
    labs(title = "Number of Variants Private to Reference Along Chromosomes",
         subtitle= paste("Number of variants : ", scales::comma(n_vars)),
         y = "Number of variants")+
    theme_classic()+
    theme(legend.position = "none")

print("Saving plot...")
ggsave(lineage_barplot_path, g, width = 16, height = 9)
ggsave(ref_private_path, p, width = 16, height = 9)

print("Done!")