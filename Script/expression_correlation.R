#!/usr/bin/env Rscript

######################################################################
# Correlation Analysis of Isoform Expression Data
######################################################################
#
# This script performs a correlation analysis on gene expression data,
# filtering based on a specified TPM threshold. It is designed for
# comparisons within the same sequencing technology
#
# The input data frame should have isoform names as rows, sample names as columns,
# and TPM values as the data.
#
# Usage:
# The script expects three command-line arguments:
# 1. Path to the input file containing the expression data.
# 2. Output file path where the correlation results will be saved.
# 3. (Optional) TPM threshold for filtering genes. Default is 0.1.
#
# Example:
# Rscript this_script.R <inputFilePath> <outputFilePath> [TPMThreshold]
#
# Author: Mahoko T. Ueda
# Date: March 29, 2024
#
######################################################################

args <- commandArgs(trailingOnly = TRUE)

# Validate minimum arguments
if(length(args) < 2) {
  stop("Not enough arguments. Usage: Rscript this_script.R <inputFilePath> <outputFilePath> [TPMThreshold]")
}

# Assign command-line arguments to variables with default values for optional arguments
inputFilePath <- args[1]
outputFilePath <- args[2]
TPMThreshold <- ifelse(length(args) >= 3, as.numeric(args[3]), 0.1)

library(readr)
library(dplyr)
library(Hmisc)

# Read data
df <- read.table(inputFilePath, sep="\t", header=TRUE) %>%
  as.data.frame()

# Apply TPM threshold filter
df_filtered <- df %>%
  filter(rowSums(. > TPMThreshold, na.rm = TRUE) > 0)

calculate_cor_pvalue <- function(filt) {
  n <- ncol(filt)
  results <- data.frame(Sample1 = character(), Sample2 = character(), Correlation = numeric(), P_Value = numeric(), R_Squared = numeric())
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      cor_test <- cor.test(filt[,i], filt[,j], method = "pearson")
      r_value <- cor_test$estimate 
      r_squared <- r_value^2
      results <- rbind(results, data.frame(Sample1 = colnames(filt)[i], Sample2 = colnames(filt)[j], Correlation = r_value, R_Squared = r_squared, P_Value = cor_test$p.value))
    }
  }
  
  return(results)
}

# Calculate correlation and p-value
results <- calculate_cor_pvalue(df_filtered)

# Save results
write.table(results, outputFilePath, sep="\t", quote=FALSE, row.names=FALSE)

cat("Calculation completed. Results saved to:", outputFilePath, "\n")

