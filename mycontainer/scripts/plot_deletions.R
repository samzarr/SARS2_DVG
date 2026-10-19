#!/usr/bin/env Rscript
# ==============================================================================
# SARS-CoV-2 Deletion Variants: Genomic Coordinate 2D Scatter Plot
# ==============================================================================
# Generates the exact Start Position vs Stop Position deletion map as presented
# in Figure 1 and Figure 7B of Jiang, Dewey, & Yin (2025):
# "Deletion detection in SARS-CoV-2 genomes using multiplex-PCR sequencing..."
# (https://www.medrxiv.org/content/10.1101/2025.04.15.25325794v2)
#
# Usage:
#   Rscript plot_deletions.R <input_file_or_dir> [output_image]
#
# Examples:
#   Rscript plot_deletions.R filtered_ERR7512327_0.01_5_20_2_34.txt
#   Rscript plot_deletions.R ./results deletion_map.png
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
input_target <- ifelse(length(args) >= 1, args[1], ".")
output_image <- ifelse(length(args) >= 2, args[2], "deletion_variant_map.png")

# 1. Discover and load input files
files <- c()
if (dir.exists(input_target)) {
  # Search for filtered pipeline summary tables or raw count tables
  files <- list.files(path = input_target, pattern = ".*_0\\.01_5_20_2_34\\.txt$", full.names = TRUE)
  if (length(files) == 0) {
    files <- list.files(path = input_target, pattern = ".*count\\.txt$", full.names = TRUE)
  }
} else if (file.exists(input_target)) {
  files <- c(input_target)
} else {
  # Try glob pattern
  files <- Sys.glob(input_target)
}

if (length(files) == 0) {
  stop(paste("Error: No deletion variant files found for target:", input_target))
}

cat(paste("Loading deletion variant data from", length(files), "file(s)...\n"))

data_list <- lapply(files, function(f) {
  df <- read_tsv(f, show_col_types = FALSE)
  if (nrow(df) > 0) {
    df$sample_id <- sub("(_annotated.*|_0\\.01.*|\\.txt$)", "", basename(f))
  }
  df
})

data_combined <- bind_rows(data_list)

if (nrow(data_combined) == 0) {
  stop("Input files contain 0 deletion records.")
}

# 2. Prepare coordinates and frequency calculations if raw count tables were provided
if (!"Frequency" %in% colnames(data_combined)) {
  if ("MinCov" %in% colnames(data_combined)) {
    data_combined$Frequency <- data_combined$depth / data_combined$MinCov
  } else {
    data_combined$Frequency <- 0.01
  }
}

if (!"logFreq" %in% colnames(data_combined)) {
  data_combined$logFreq <- log10(data_combined$Frequency)
}

data_combined$logFreq[!is.finite(data_combined$logFreq)] <- -2
data_forward <- data_combined %>% filter(start < end)

cat(paste("Plotting", nrow(data_forward), "deletion variant junction(s)...\n"))

# 3. Build Figure 1 reproduction plot using author's original ggplot parameters
p <- ggplot(data_forward, aes(x = end, y = start, alpha = logFreq)) +
  # 5' to 3' diagonal boundary reference (y = x)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray70", size = 0.5) +
  # Deletion junction points (Stop on X, Start on Y)
  geom_point(size = 1.2, color = "blue") +
  theme_linedraw(base_size = 16) +
  xlim(0, 31500) +
  ylim(0, 31500) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.position = c(0.15, 0.65),
    legend.background = element_rect(fill = "white", color = "gray80"),
    aspect.ratio = 1
  ) +
  labs(
    x = "Stop Position (nt)",
    y = "Start Position (nt)",
    title = "SARS-CoV-2 Deletion Variants",
    subtitle = paste("Mapped junctions across", length(unique(data_forward$sample_id)), "sample(s)")
  ) +
  scale_alpha_continuous(
    limits = c(-2, 0),
    breaks = seq(-2, 0, by = 0.5),
    labels = c("-2.0", "-1.5", "-1.0", "-0.5", "0")
  )

# 4. Save output
ggsave(
  filename = output_image,
  plot = p,
  width = 6,
  height = 6,
  dpi = 300
)

cat(paste("Successfully generated deletion variant map:", output_image, "\n"))
