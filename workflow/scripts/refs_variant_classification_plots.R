log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

print("Loading libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(patchwork))

print("Reading input parameters...")

chromosomes_path = snakemake@input$chromosomes
variant_classification_path = snakemake@input$variants
presence_path = snakemake@input$presence

lineage=snakemake@wildcards$lineage
window_size = snakemake@params$window_size

lineage_barplot_path = snakemake@output$plot
ref_private_path = snakemake@output$plot_density

print("Reading files...")
vars_classification = read.delim(variant_classification_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
presence = read.delim(presence_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
chromosomes = read.delim(chromosomes_path, sep = ",", header = TRUE, stringsAsFactors = TRUE)

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
p <- ggplot(status_density)+
    geom_segment(aes(x=1, xend=length, y=-1, yend=-1), linewidth = 0.1, color = "black")+
    geom_col(aes(x = window_start, y = n), color = "black", fill = "black")+
    facet_wrap(~chromosome, ncol = 2, strip.position = "right")+
    scale_x_continuous(name = "Position (bp) ", labels = comma)+
    labs(title = "Number of Variants Private to Reference Along Chromosomes",
        subtitle= paste("Number of variants : ", scales::comma(n_vars)),
            y = "Number of variants")+
    theme_classic()+
    theme(legend.position = "none")


print("Saving plot...")
ggsave(lineage_barplot_path, c, width = 8, height = 5) # Choose g, b, or c
ggsave(ref_private_path, p, width = 16, height = 9)

print("Done!")