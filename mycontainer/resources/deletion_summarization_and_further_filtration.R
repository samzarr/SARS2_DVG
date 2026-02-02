#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(tibble)

# -----------------------------
# STEP 1: Summarize recombination data
# -----------------------------
recomb_files <- list.files(pattern = "*deletion.sorted.txt")

process_file <- function(file) {
  data <- read_tsv(file)
  
  summarized_data <- data %>%
  group_by(chrom, start, end) %>%
  summarize(depth = n(), 
        depth_positive = sum(strand == '+'),
        depth_negative = sum(strand == '-'), 
        n_fragments = n_distinct(read_name), 
        max_left_overhang = max(left_overhang), 
        max_right_overhang = max(right_overhang), 
        max_aligned_length = max(aligned_length)) %>%
  ungroup()
  
  output_file <- sub("deletion.sorted.txt$", "count.txt", file)
  write_tsv(summarized_data, file = output_file)
}

lapply(recomb_files, process_file)

# -----------------------------
# STEP 2: Frequency calculation and filtering
# -----------------------------
coverage_files <- list.files(path = ".", pattern = "*out.sorted.coverage", full.names = TRUE)
count_files <- list.files(path = ".", pattern = "*count.txt", full.names = TRUE)

for (bed_file in count_files) {
  base_name <- sub("_annotated_standardizationAligned\\.count\\.txt$", "", basename(bed_file))
  coverage_file <- coverage_files[grepl(base_name, coverage_files)]
  
  if (length(coverage_file) == 0) {
    print(paste("No coverage file found for:", base_name))
    next
  }
  
  data <- read.table(bed_file, sep = "", header = TRUE)
  data <- data[order(data$depth), ]
  data <- data[which(data$end < 29903), ]
  junction_depth <- sum(data$depth)
  print(paste0("The total junction depth is: ", junction_depth))
  
  coverage <- read.table(coverage_file, sep = "", header = FALSE)
  coverage <- coverage %>% rename(Genome = V1, Position = V2, Depth = V3)
  total_depth <- sum(coverage$Depth)
  print(paste0("The total mapped depth is: ", total_depth))
  
  data_forward <- data[which(data$start < data$end), ]
  data_forward$start_cov <- coverage[unlist(data_forward[2]), 3]
  data_forward$end_cov <- coverage[unlist(data_forward[3] + 1), 3]
  data_forward <- transform(data_forward, MinCov = pmin(start_cov, end_cov))
  data_forward$Frequency <- data_forward$depth / data_forward$MinCov
  data_forward$logFreq <- log10(data_forward$Frequency)
  data_forward$logFreq[!is.finite(data_forward$logFreq)] <- -2
  data_forward$deletion_length <- data_forward$end - data_forward$start
  
  
  data_forward <- data_forward[
    data_forward$depth >= 5 &
    data_forward$Frequency >= 0.01 &
    data_forward$MinCov > 20 &
    data_forward$depth_positive > 2 &
    data_forward$depth_negative > 2 &
    data_forward$max_right_overhang > 34 &
    data_forward$max_left_overhang > 34,
  ]

  #data_forward <- data_forward[
    #data_forward$depth >= 1 &          
    #data_forward$Frequency >= 0 &  
    #data_forward$MinCov > 0 &
    #data_forward$max_right_overhang >= 0 &
    #data_forward$max_left_overhang >= 0,
  #]
  
  cat("Checking unique deletion coordinates:\n")
  print(unique(data_forward[, c("start", "end")]))
  cat("Total records:", nrow(data_forward), "\n")
  cat("Unique (start, end) pairs:", nrow(unique(data_forward[, c('start','end')])), "\n\n")


  write.table(data_forward, file = paste0(base_name, "_0.01_5_20_2_34.txt"), sep = "\t", row.names = FALSE)
}
