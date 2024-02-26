log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))

regions <-  read.table(snakemake@input[[1]], header = TRUE, stringsAsFactors = TRUE, sep = "\t")

size_threshold <- 1000
change_threshold <- 0.5
size_threshold <- snakemake@params[[1]]
change_threshold <- snakemake@params[[2]]
region_size <- regions$End[1]

dupdels <- list()
for (chrom in levels(as.factor(regions$Accession))){
  regions_select <- regions %>%
    select(Accession, Start, End, Smooth, Norm_Median)%>%
    filter(Accession == chrom )

  regions_windowed <- regions_select
  regions_windowed$Structure <- "Haploid"
  for (i in 1:nrow(regions_windowed)){
    if(regions_windowed$Smooth[i] >= 1 + change_threshold){
      regions_windowed$Structure[i] <- "Duplication"
    } else if (regions_windowed$Smooth[i] <= 1 - change_threshold) {
      regions_windowed$Structure[i]<- "Deletion"
    } else {
      regions_windowed$Structure[i] <- "Haploid"
    }
  }
  
  regions_windowed$window_index <- 1
  for (i in 2:nrow(regions_windowed)){
    regions_windowed$window_index[i] <- ifelse(regions_windowed$Structure[i]  == regions_windowed$Structure[i-1],
                                          regions_windowed$window_index[i-1],
                                          regions_windowed$window_index[i-1]+1 )
  }
  
  windows<- regions_windowed %>%
    group_by(window_index)%>%
    mutate(Window_Smooth_Cov = round(mean(Smooth), 2),Window_Norm_Cov = round(mean(Norm_Median)), n = n())%>%
    mutate(Win_Start = Start[1], Win_End = End[n])%>%
    filter(row_number()==1)%>%
    ungroup()%>%
    mutate(Window_Size = n*region_size)%>%
    select(-c(Smooth, window_index, Start, End))%>%
    select(Accession, Start = Win_Start, End = Win_End, Window_Size, Window_Norm_Cov, Window_Smooth_Cov, Structure)
  

  chrom_dupdels <- windows %>%
    filter(Window_Size >= size_threshold)%>%
    filter(Structure != "Haploid")

  dupdels[[chrom]] <- chrom_dupdels
}

dupdels_df <- do.call(rbind,dupdels)

write.table(dupdels_df, snakemake@output[[1]], col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
