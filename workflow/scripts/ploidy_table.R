log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))

regions <-  read.csv(snakemake@input[[1]], header = TRUE, stringsAsFactors = TRUE)

diff_threshold <- snakemake@params[[1]]
size_threshold <- snakemake@params[[2]]
change_threshold <- snakemake@params[[3]]
region_size <- regions$End[1]

dupdels <- list()
for (chrom in levels(as.factor(regions$Chromosome))){
  regions_select <- regions %>%
    select(Chromosome, Start, End, Norm_Median)%>%
    filter(Chromosome == chrom)
  

  regions_windowed <- regions_select
  regions_windowed$window_index <- 1
  for (i in 2:nrow(regions_windowed)){
    regions_windowed$window_index[i] <- ifelse(abs(regions_windowed$Norm_Median[i]-regions_windowed$Norm_Median[i-1]) < diff_threshold,
                                          regions_windowed$window_index[i-1],
                                          regions_windowed$window_index[i-1]+1 )
  }
  
  windows<- regions_windowed %>%
    group_by(window_index)%>%
    mutate(window_mean = round(mean(Norm_Median), 2), n = n())%>%
    filter(row_number()==1)%>%
    ungroup()%>%
    select(-c(Norm_Median, window_index))
  
  windows$End <- 0
  for (i in 1:nrow(windows)){
    windows$End[i]<- windows$Start[i+1]
  } 
  windows$Structure <- "Haploid"
  for (i in 1:nrow(windows)){
    if(windows$window_mean[i] > 1 + change_threshold){
       windows$Structure[i] <- "Duplication"
    } else if (windows$window_mean[i] < 1 - change_threshold) {
      windows$Structure[i]<- "Deletion"
    } else {
      windows$Structure[i] <- "Haploid"
    }
  }
  
  chrom_dupdels <- windows %>%
    mutate(window_size = n*region_size)%>%
    filter(window_size >= size_threshold)%>%
    filter(Structure != "Haploid")%>%
    select(-n)
  
  dupdels[[chrom]] <- chrom_dupdels
}

dupdels_df <- do.call(rbind,dupdels)

write_csv(dupdels_df, snakemake@output[[1]], col_names = TRUE)
