#!/usr/bin/env Rscript

######################################################################
# Comparative Analysis of Gene Expression Between Short-read and Long-read Data
######################################################################
#
# This script compares gene expression data obtained from short-read and long-read
# sequencing technologies. It filters genes based on a TPM threshold and calculates
# correlation between expression levels across the technologies.
#
# Usage:
# The script expects the following command-line arguments:
# 1. Path to the long-read expression data file.
# 2. Path to the short-read expression data file.
# 3. Output directory path for saving filtered data and results.
# 4. (Optional) TPM threshold for filtering genes. Default is 0.1.
#
# Example:
# Rscript this_script.R <longReadFilePath> <shortReadFilePath> <outputDirPath> [TPMThreshold]
#
# Dependencies:
# - readr, dplyr, ggplot2 R packages
#
# It is designed for comparisons within the same set of samples across different 
# sequencing technologies, focusing on data consistency.
#
# Author: Mahoko T. Ueda
# Date: March 29, 2024
#
######################################################################

args <- commandArgs(trailingOnly = TRUE)

# Validate minimum arguments
if(length(args) < 3) {
  stop("Not enough arguments. Usage: Rscript this_script.R <longReadFilePath> <shortReadFilePath> <outputDirPath> [TPMThreshold]")
}

# Assign command-line arguments to variables with default values for optional arguments
longReadFilePath <- args[1]
shortReadFilePath <- args[2]
outputDirPath <- args[3]
TPMThreshold <- ifelse(length(args) >= 4, as.numeric(args[4]), 0.1)

library(readr)
library(dplyr)
library(ggplot2)

# Read data
long_read_df <- read.table(longReadFilePath, sep="\t", header=TRUE)
short_read_df <- read.table(shortReadFilePath, header=TRUE, sep="\t")

# Select samples present in both datasets
selected_samples <- intersect(colnames(long_read_df), colnames(short_read_df))
long_read_sel <- long_read_df %>% select(all_of(selected_samples))
short_read_sel <- short_read_df %>% select(all_of(selected_samples))

# Extract shared genes
common_genes_names <- intersect(rownames(long_read_sel), rownames(short_read_sel))
long_common <- long_read_sel[common_genes_names, ]
short_common <- short_read_sel[common_genes_names, ]

# Filtering: keep genes with TPM > threshold in all samples for both technologies
filter_long <- apply(long_common, 1, function(x) all(x > TPMThreshold))
filter_short <- apply(short_common, 1, function(x) all(x > TPMThreshold))
final_common_genes <- filter_long & filter_short

# Filtered data
long_filtered <- long_common[final_common_genes, ]
short_filtered <- short_common[final_common_genes, ]

# Write filtered data to files
write.table(long_filtered, file.path(outputDirPath, "long-read_filtered_TPM0.1.txt"), sep="\t", quote=FALSE)
write.table(short_filtered, file.path(outputDirPath, "short-read_filtered_TPM0.1.txt"), sep="\t", quote=FALSE)

# Log conversion
long_log <- log2(long_filtered + 1)
short_log <- log2(short_filtered + 1)

# Calculate correlation for each sample
samples <- colnames(long_filtered)
cor_results <- sapply(samples, function(sample) cor(long_log[[sample]], short_log[[sample]], use = "complete.obs"))

# Display sample names and correlation coefficients
names(cor_results) <- samples
print(cor_results)

# Calculate and display mean and standard deviation of correlations
mean_corr <- mean(cor_results, na.rm = TRUE)
std_corr <- sd(cor_results, na.rm = TRUE)
print(paste("Mean correlation:", mean_corr))
print(paste("Standard deviation:", std_corr))
