##########################################################################
# Simulating RNA-seq Data with Polyester Based on Quantification Results
##########################################################################
#
# This script simulates RNA-seq data based on Salmon quantification results,
# using the Polyester package. 
#
# Usage:
# The script expects three to five command-line arguments:
# 1. Path to the input directory containing Salmon's quant.sf files.
# 2. Output directory path for storing simulated data.
# 3. Path to the transcriptome fasta file used for simulation.
# 4. (Optional) Number of cores for parallel processing. Default is 5.
# 5. (Optional) Read length for simulation. Default is 150.
# 6. (Optional) Error rate for simulation. Default is 0.00109.
#
# Example:
# Rscript simulate_rnaseq.R <inputDirPath> <outputDirPath> <fastaPath> [numCores] [readLength] [errorRate]
#
# Author: Mahoko T. Ueda
# Date: March 27, 2024
#
##########################################################################

args <- commandArgs(trailingOnly = TRUE)

# Validate arguments
if(length(args) != 3) {
  stop("Not enough arguments. Usage: Rscript simulate_rnaseq.R <inputDirPath> <outputDirPath> <timePointName>")
}

# Assign command-line arguments to variables
inputDirPath <- args[1]
outputDirPath <- args[2]
fasta_path <- args[3]
mc.cores <- ifelse(length(args) >= 4, as.numeric(args[4]), 5)
read_length <- ifelse(length(args) >= 5, as.numeric(args[5]), 150)
error_rate <- ifelse(length(args) >= 6, as.numeric(args[6]), 0.00109)


# Error rate for NovaSeq; adjust based on your data
error_rate <- 0.00109

library(polyester)
library(readr)
library(Biostrings)
library(parallel)

# Ensure output directory exists
if (!dir.exists(outputDirPath)) {
  dir.create(outputDirPath, recursive = TRUE)
}

# Obtain list of sample directories
sample_names <- list.files(path = inputDirPath)

# Process each sample
results <- mclapply(sample_names, function(sample_name) {
  print(paste("Processing sample:", sample_name))
  
  quant_sf_path <- file.path(inputDirPath, sample_name, "quant.sf")
  quant_data <- read_tsv(quant_sf_path, col_names = TRUE, show_col_types = FALSE)
  
  expression_levels <- round(quant_data$NumReads)
  expression_levels[expression_levels == 0] <- 1
  num_transcripts <- length(expression_levels)
  fold_changes <- matrix(rep(1, num_transcripts), nrow = num_transcripts, ncol = 1)
  
  output_dir <- file.path(outputDirPath, sample_name)
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  simulate_experiment(
    fasta = fasta_path,
    error_model = "uniform",
    error_rate = error_rate,
    reads_per_transcript = expression_levels,
    outdir = output_dir,
    fold_changes = fold_changes,
    paired = TRUE,
    num_reps = c(1),
    readlen = 150
  )
}, mc.cores = mc.cores)  # Adjust for your system

print("Simulation completed.")
