#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)

# -----------------------------
# STEP 1: Summarize recombination data
# -----------------------------
if (length(args) >= 1 && file.exists(args[1])) {
  recomb_files <- args[1]
} else {
  recomb_files <- list.files(pattern = ".*deletion\\.sorted\\.txt$")
}

process_file <- function(file) {
  data <- read_tsv(file, show_col_types = FALSE)
  
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
  
  output_file <- sub("deletion\\.sorted\\.txt$", "count.txt", file)
  write_tsv(summarized_data, file = output_file)
}

lapply(recomb_files, process_file)

# -----------------------------
# STEP 2: Frequency calculation and filtering
# -----------------------------
coverage_files <- list.files(path = ".", pattern = ".*out\\.sorted\\.coverage$", full.names = TRUE)
count_files <- list.files(path = ".", pattern = ".*count\\.txt$", full.names = TRUE)

if (length(args) >= 1 && !file.exists(args[1])) {
  # If a sample name prefix was passed as an argument
  count_files <- count_files[grepl(args[1], count_files)]
}

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
  
  coverage <- read.table(coverage_file[1], sep = "", header = FALSE)
  coverage <- coverage %>% rename(Genome = V1, Position = V2, Depth = V3)
  total_depth <- sum(coverage$Depth)
  print(paste0("The total mapped depth is: ", total_depth))
  
  data_forward <- data[which(data$start < data$end), ]
  start_idx <- match(data_forward$start, coverage$Position)
  end_idx <- match(data_forward$end + 1, coverage$Position)
  data_forward$start_cov <- ifelse(is.na(start_idx), 0, coverage$Depth[start_idx])
  data_forward$end_cov <- ifelse(is.na(end_idx), 0, coverage$Depth[end_idx])
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
  
  plot <- ggplot(data_forward, aes(end, start, alpha = logFreq)) + 
    geom_point(size = 1, color = "blue") + 
    theme_linedraw(base_size = 20) + 
    ylim(0, 31500) + xlim(0, 31500) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 11),
          legend.position = c(0.15, 0.63)) + 
    labs(x = "Stop Position (nt)", y = "Start Position (nt)") +
    scale_alpha_continuous(limits = c(-2, 0), 
                           breaks = seq(-2, 0, by = 0.5), 
                           labels = c("-2.0", "-1.5", "-1.0", "-0.5", "0"))
  
  print(plot)
  ggsave(filename = paste0(base_name, "_0.01_5_20_2_34_blue.tiff"), plot = plot,
         scale = 1, width = 4.5, height = 4, units = "in", dpi = 1200, limitsize = TRUE)
}
