#!/usr/bin/env /usr/local/apps/R/4.5/4.5.2/bin/Rscript

library(tidyverse)

# Get command-line arguments (trailingOnly=TRUE ignores Rscript arguments)
args <- commandArgs(trailingOnly = TRUE)

input_file <- if(is.na(args[1])) {
  stop("Please provide an input file as the first argument.")
} else {
  args[1]
}
output_file <- if(is.na(args[2])) {
  "filtered_positions.bed"
} else {
  str_replace(args[2], pattern = "\\.txt", replacement = ".bed")
}

# Read dataset
#dp_data <- read_tsv("data/genome/high_qual_snps_AD5_AF0.9.filtered.DP.txt",
dp_data <- read_tsv(args[1], col_names = c("CHROM", "POS"))

# Convert to data.table to improve performance
DT <- dp_data %>%
  pivot_longer(cols = -c(CHROM, POS), names_to = "tmp", values_to = "DP") %>%
  data.table::as.data.table()

DT[, c("strain", "DP") := data.table::tstrsplit(DP, split = ":", fixed = TRUE)]

# Set DP to integer
DT[, DP := as.integer(DP)]

# Calculate Mean and SD by strain
DT[, c("Mean", "SD") := .(mean(DP, na.rm = TRUE), sd(DP, na.rm = TRUE)), by = strain] 

# Set CNV to true if DP is greater than mean + 2*SD
DT[, CNV := DP > (Mean + 2 * SD)]

# Filter out genomic positions where CNV == TRUE in any strain
DT[, Filtered := any(CNV), by = .(CHROM, POS)]

# Add start column for BED format (0-based)
DT[, Start := POS - 1]

# Keep only rows where Filtered is FALSE
DT[Filtered == FALSE, .(CHROM, Start, POS)] %>%
  distinct() %>%
  write_tsv(output_file, col_names = FALSE)