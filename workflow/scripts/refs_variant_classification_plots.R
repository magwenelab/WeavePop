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

chromosomes_path <- snakemake@input$chromosomes
variant_classification_path <- snakemake@input$classif
variants_path <- snakemake@input$variants
presence_path <- snakemake@input$presence
loci_path <- snakemake@input$loci

lin<-snakemake@wildcards$lineage
window_size <-snakemake@params$window_size

lineage_vars_path <- snakemake@output$plot
ref_private_path <- snakemake@output$plot_density
   

print("Reading files...")
vars_classification <- read.delim(variant_classification_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
variants <- read.delim(variants_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
presence <- read.delim(presence_path, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
chromosomes <- read.delim(chromosomes_path, sep = ",", header = TRUE, stringsAsFactors = TRUE)
loci_table <- read.delim(loci_path, header = TRUE, stringsAsFactors = TRUE, na = c("", "N/A"))

vars_classification$category <- factor(vars_classification$category, levels = c("Reference private","Private", "Non-private"))
vars_classification$impact<- factor(vars_classification$impact, levels = c("High","Moderate","Low", "Modifier"))

vars_classification <- left_join(variants, vars_classification, by = "var_id")

print("Ordering chromosome names...")
chromosomes <- chromosomes %>%
    filter(lineage == lin)

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
    select(var_id, lineage, chromosome, impact, category )

total_vars = nrow(vars_type)

vars_type_sample <- left_join(presence, vars_type, by = "var_id")%>%
    group_by(sample, lineage, category, impact)%>%
    summarize(n_variants = n())%>%
    ungroup()

vars_type_sample_max <- vars_type_sample %>%
    group_by(category)%>%
    summarize(max_variants = max(n_variants) + max(n_variants)/10)

vars_summary <- vars_type %>%
    group_by(impact, category)%>%
    summarize(n_variants = n())%>%
    left_join(vars_type_sample_max, by = "category")

g <- ggplot(vars_type_sample, aes(x = impact, y = n_variants))+
    expand_limits(y = 0)+
    geom_quasirandom(aes(color = impact), alpha = 0.5)+
    geom_boxplot(alpha = 0)+
    geom_text(data = vars_summary, aes(y = max_variants, x = impact, label = scales::comma(n_variants)))+
    facet_wrap(~category, nrow = 1, scales = "free")+
    scale_y_continuous(name = "Number of variants per sample", labels = comma)+
    labs(title = paste("Number of Variants per Impact and Privateness Category in each Sample of Lineage ", lin),
         subtitle = paste("Total variants:", scales::comma(total_vars)),
         x = "", 
         color = "Impact")+
    theme_classic()
print("Saving plot...")
ggsave(lineage_vars_path, g, width = 16, height = 9)

print("Making density plot of reference private variants")

if (filter(vars_classification, category == "Reference private") %>% nrow() == 0) {
  print("No reference private variants, creating empty plot...")
  u <- ggplot() + 
        theme_classic() +
        labs(title = paste("No reference private variants found in", lin))
  ggsave(ref_private_path, u, width = 4, height = 1)
  print("Done!")
  quit(status = 0)
}

vars_windows <- vars_classification %>%
    left_join(chromosomes, by = c("accession", "lineage"))%>%
    group_by(chromosome) %>%
    mutate(window = floor(pos / window_size) + 1)%>%
    ungroup()%>%
    group_by(chromosome, window) %>%
    mutate(window_start = (window - 1) * window_size + 1,
            window_end = window * window_size)%>%
    ungroup()%>%
    filter(category == "Reference private")

category_density <- vars_windows %>% 
    group_by(chromosome, window,window_start, window_end, category, length)%>%
    summarize(n = as.integer(n()))

n_vars <- length(unique(vars_windows$var_id))

print("Plotting category density")
p <- ggplot()+
    geom_segment(data = chromosomes, 
                aes(x=1, xend=length, y=0, yend=0), 
                linewidth = 1, color = "black", alpha = 0.1)+
    geom_col(data = category_density, 
            aes(x = window_start, y = n), 
            color = "black", fill = "black")+
    facet_wrap(~chromosome, ncol = 2, strip.position = "right")+
    scale_x_continuous(name = "Position (bp) ", labels = comma)+
    scale_y_continuous(name = "Number of variants", # Use only integers
                       breaks = function(limits) {
                         b <- pretty(limits)
                         b <- round(b)
                         b <- unique(b)
                         b <- b[b >= 0]
                         b
                       },
                       labels = comma)+
    labs(title = "Number of Variants in Windows Along Chromosomes Private to the Reference Genome",
         subtitle= paste("Lineage: ", lin, " Window size: ", window_size, " Total variants: ", scales::comma(n_vars)))+
    theme_classic()

max_vars <- max(category_density$n)

print("Adding loci data if available...")
if (nrow(loci_table)!= 0){
  loci_sample <- loci_table %>% 
    filter(lineage == lin) %>%
    select(accession, start , end , loci)%>%
    droplevels()
  loci <- left_join(loci_sample, chromosomes, by = c("accession"))%>%
    select(accession, chromosome, start, end, loci)
  dark2 <- brewer.pal(8, "Dark2")[1:6]
  l_colors <- dark2[1:nlevels(loci$loci)]
  l_lim <- max_vars + max_vars*0.1
  print("Adding loci to plot...")
  p <- p +  
      geom_point(data = loci, aes(x=start, y = l_lim, color = loci))+  
      scale_color_manual(name = "Features", values = l_colors)+
      guides(color = guide_legend(order=3))
}

print("Saving plot...")
ggsave(ref_private_path, p, width = 16, height = 9)

print("Done!")