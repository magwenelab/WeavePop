library(tidyverse)
library(scales)
library(ggbeeswarm)
library(patchwork)

effects_path  = "/FastData/czirion/WeavePop/test_bigger/results/04.Intermediate_files/02.Dataset/snpeff/VNI_effects.tsv"
presence_path = "/FastData/czirion/WeavePop/test_bigger/results/04.Intermediate_files/02.Dataset/snpeff/VNI_presence.tsv"
variants_path = "/FastData/czirion/WeavePop/test_bigger/results/04.Intermediate_files/02.Dataset/snpeff/VNI_variants.tsv"
metadata_path = "/FastData/czirion/WeavePop/test_bigger/results/02.Dataset/metadata.csv"
chromosomes_path = "/FastData/czirion/WeavePop/test_bigger/results/02.Dataset/chromosomes.csv"
lineage = "VNBI"
lineage_barplot_path = paste("./", lineage, "_variants_summary_barplot.png", sep = "")
plots_path = "/FastData/czirion/WeavePop/test_bigger/results/01.Samples/plots/"

variants = read.delim(variants_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
effects = read.delim(effects_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
presence = read.delim(presence_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
metadata = read.delim(metadata_path, sep = ",", header = TRUE, stringsAsFactors = TRUE)
chromosomes = read.delim(chromosomes_path, sep = ",", header = TRUE, stringsAsFactors = TRUE)

# Add window information to variants table
n <- 1000

vars_windows <- variants %>%
    left_join(chromosomes, by = c("accession", "lineage"))%>%
    group_by(chromosome) %>%
    mutate(window = floor(pos / n) + 1)%>%
    ungroup()%>%
    group_by(chromosome, window) %>%
    mutate(window_start = (window - 1) * n + 1,
            window_end = window * n)%>%
    ungroup()
    
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

vars_windows_effs_status <- left_join(vars_windows, effs, by = "var_id") %>%
    left_join(samples_in_lineage, by = "lineage") %>%
    left_join(presence_summary, by = "var_id") %>%
    mutate(status = ifelse(n_samples == samples_in_lineage, "Reference private", ifelse(n_samples == 1, "Private", "Non-private") ))%>%
    select(-c(samples_in_lineage, n_samples))
vars_windows_effs_status$status <- factor(vars_windows_effs_status$status, levels = c("Reference private","Private", "Non-private"))

# # Average number of private variants per sample
# vars_windows_effs_status %>%
#     group_by(lineage,status)%>%
#     summarize(n_variants = length(unique(var_id)))%>%
#     filter(status == "Private")%>%
#     left_join(samples_in_lineage, by = "lineage")%>%
#     mutate(avg_per_sample = n_variants / samples_in_lineage)%>%
#     pull(avg_per_sample)

summary_effs_status <- vars_windows_effs_status %>%
    group_by(status, impact)%>%
    summarize(n_variants = length(unique(var_id)))

vars_type <- vars_windows_effs_status%>%
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

# Barplot per sample

for (samp in levels(presence$sample)){

    vars_sample <- presence%>%
        filter(sample == samp)%>%
        droplevels()
    
    n_vars <- length(unique(vars_sample$var_id))

    vars <- vars_windows_effs_status%>%
        filter(var_id %in% vars_sample$var_id,
                status != "Reference private")
    
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

    path_e <- paste(plots_path, samp, "_impacts.png", sep = "")
    path_p <- paste(plots_path, samp, "_status.png", sep = "")
    path_b <- paste(plots_path, samp, "_bar.png", sep = "")

    ggsave(path_e, e, width = 16, height = 9)
    ggsave(path_p, p, width = 16, height = 9)
    ggsave(path_b, b, width = 8, height = 5)

}



# Get number of variants with each type of impact per window

# effs_density <- vars_windows_effs_status%>% 
#     group_by(chromosome, window,window_start, window_end)%>%
#     summarize(Modifier = sum(MODIFIER),
#                 Low = sum(LOW),
#                 Moderate = sum(MODERATE),
#                 High = sum(HIGH))%>%
#     pivot_longer(cols = c("High", "Moderate","Low","Modifier"),
#                 names_to = "impact",
#                 values_to = "n_variants")

# e <- ggplot(effs_density)+
#     geom_col(aes(x = window_start, y = n_variants, color = impact))+
#     facet_grid(chromosome~impact)+#, ncol =1, strip.position = "right")+
#         labs(title="Number of variants per window with at least one effect per impact type",
#                 y = "Number of variants")+
#     theme(panel.grid = element_blank(),
#             panel.grid.major.x = element_blank(),
#             panel.grid.minor.x = element_blank(),
#             panel.grid.major.y = element_blank(),
#             panel.grid.minor.y = element_blank(),
#             panel.background = element_blank(),
#             panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 2),
#             legend.position = "none")

# e

ggsave(lineage_barplot_path, b, width = 8, height = 5) # Choose g, b, or c
