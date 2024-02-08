log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))

regions <-  read.csv(snakemake@input[[1]], header = FALSE,sep = "\t", col.names = c("Chromosome", "Start", "Coverage", "Smooth"), colClasses = c("factor", "numeric", "numeric", "numeric"))
regions <- regions %>%
  mutate(End = Start + 500)
head(regions)

diff_threshold <- snakemake@params[[1]]
size_threshold <- snakemake@params[[2]]
change_threshold <- snakemake@params[[3]]
region_size <- regions$End[1]

dupdels <- list()
for (chrom in levels(as.factor(regions$Chromosome))){
  regions_select <- regions %>%
    filter(Chromosome == chrom) 
  

  regions_windowed <- regions_select
  regions_windowed$window_index <- 1
  for (i in 2:nrow(regions_windowed)){
    regions_windowed$window_index[i] <- ifelse(abs(regions_windowed$Smooth[i]-regions_windowed$Smooth[i-1]) < diff_threshold,
                                          regions_windowed$window_index[i-1],
                                          regions_windowed$window_index[i-1]+1 )
  }
  
  print(regions_windowed)

  windows<- regions_windowed %>%
    group_by(window_index)%>%
    mutate(Window_Norm_Cov = round(mean(Smooth), 2), n = n())%>%
    mutate(Win_Start = Start[1], Win_End = End[n])%>%
    filter(row_number()==1)%>%
    ungroup()%>%
    mutate(Window_Size = n*region_size)%>%
    select(-c(Smooth, window_index, Start, End))%>%
    select(Chromosome, Start = Win_Start, End = Win_End, Window_Size, Window_Norm_Cov )
  
  windows$Structure <- "Haploid"
  for (i in 1:nrow(windows)){
    if(windows$Window_Norm_Cov[i] > 1 + change_threshold){
       windows$Structure[i] <- "Duplication"
    } else if (windows$Window_Norm_Cov[i] < 1 - change_threshold) {
      windows$Structure[i]<- "Deletion"
    } else {
      windows$Structure[i] <- "Haploid"
    }
  }
  
  chrom_dupdels <- windows %>%
    filter(Window_Size >= size_threshold)%>%
    filter(Structure != "Haploid")

  dupdels[[chrom]] <- chrom_dupdels
}

dupdels_df <- do.call(rbind,dupdels)

write_csv(dupdels_df, snakemake@output[[1]], col_names = TRUE)
