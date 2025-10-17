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
variants_path = snakemake@input$variants
presence_path = snakemake@input$presence
chromosomes_path = snakemake@input$chromosomes
barplot_path = snakemake@output$barplot
status_path = snakemake@output$status
impact_path = snakemake@output$impact

# Add window information to variants table
n <- 1000
variants = read.delim(variants_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
presence = read.delim(presence_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
chromosomes = read.delim(chromosomes_path, sep = ",", header = TRUE, stringsAsFactors = TRUE)


vars_windows <- variants %>%
    left_join(chromosomes, by = c("accession", "lineage"))%>%
    group_by(chromosome) %>%
    mutate(window = floor(pos / n) + 1)%>%
    ungroup()%>%
    group_by(chromosome, window) %>%
    mutate(window_start = (window - 1) * n + 1,
            window_end = window * n)%>%
    ungroup()
    
vars_sample <- presence %>%
    filter(sample == samp)%>%
    droplevels()

n_vars <- length(unique(vars_sample$var_id))

vars <- vars_windows %>%
    filter(var_id %in% vars_sample$var_id)

if (length(unique(presence$sample)) > 1){
    print("Removing variants private to reference genome")
    vars <- vars%>%
        filter(status != "Reference private")%>%
        droplevels()
} else {
    print("Only one sample in lineage, plotting all variants")
}

status_density <- vars %>% 
    group_by(chromosome, window,window_start, window_end, status)%>%
    summarize(n = n())

effs_density <- vars%>% 
    group_by(chromosome, window,window_start, window_end, impact)%>%
    summarize(n = n())

p <- ggplot(status_density)+
    geom_col(aes(x = window_start, y = n, color = status))+
    facet_grid(chromosome~status)+
    scale_x_continuous(name = "Position (bp) ", labels = comma)+
    labs(title = "Number of variants of each status per window ",
        subtitle= paste("Sample: ", samp," Number of variants: ", scales::comma(n_vars)),
            y = "Number of variants")+
    theme(panel.grid = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 2),
    legend.position = "none")

e <- ggplot(effs_density)+
    geom_col(aes(x = window_start, y = n, color = impact))+
    facet_grid(chromosome~impact, scales = "free_x")+#, ncol =1, strip.position = "right")+
    scale_x_continuous(name = "Position (bp) ", labels = comma)+
    labs(title="Number of variants per window per impact",
        subtitle= paste("Sample: ", samp," Number of variants: ", scales::comma(n_vars)),
            y = "Number of variants")+
    theme(panel.grid = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 2),
    legend.position = "none")

total_vars = nrow(vars)
b <- ggplot(vars, aes(x = status, fill = impact))+
    geom_bar(stat = "count", position = "dodge")+
    geom_text(stat= "count", position = position_dodge(width =0.9), 
                aes(label = scales::comma(after_stat(count))),
                vjust = -0.5,
                size = 3)+
    scale_y_continuous(name = "Number of variants", labels = comma)+
    labs(title = "Number of variants in each status and impact",
        subtitle= paste("Sample: ", samp," Total number of variants:", scales::comma(total_vars)),
        x = "", 
        fill = "Impact")+
    theme_classic()

print("Saving plots...")
ggsave(barplot_path, b, width = 16, height = 9)
ggsave(status_path, p, width = 16, height = 9)
ggsave(impact_path, e, width = 8, height = 5)
print("Done!")