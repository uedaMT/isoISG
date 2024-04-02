###############################################################
# Correlation between IR-isoform and Spliceosome Expression
###############################################################
#
# This R script performs analysis to investigate the correlation between IR-isoform and spliceosome expression.
# It calculates Pearson correlation coefficients for gene pairs, adjusts p-values for multiple testing using FDR,
# and identifies statistically significant gene pairs. Specifically, it filters out pairs where both genes are
# associated with IR events, identified by ".IR" in their names.
#
# Please ensure that you have the necessary dependencies installed before running this script.
# Author: Mahoko T. Ueda
# Date: May 30, 2023
#
# Usage: Rscript IR-spliceosome_cor-test_sc.R --input "path/to/your/input/expression_data.txt" --output "/path/to/your/output/significant_correlations_noIRpairs_FDR0.05.txt"
#
#
#----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

# Validate arguments
if(length(args) < 2) {
  stop("Not enough arguments. Usage: Rscript IR-spliceosome_cor-test_sc.R <input_file> <output_file>")
}

input_file <- args[1]
output_file <- args[2]

# Read expression data
exprs <- read.table(input_file, header=TRUE, sep="\t")
rownames(exprs) <- exprs$ID
exprs <- exprs[,-1]

# Calculate correlation coefficients
cor_matrix <- cor(t(exprs), method = "pearson")


# Function to calculate p-values
r2p <- function(r, n) {
  t_value <- r * sqrt((n - 2) / (1 - r^2))
  p_value <- 2 * pt(-abs(t_value), df = n - 2)
  return(p_value)
}

# Sample size
n <- ncol(exprs)

# Calculate p-values
p_matrix <- matrix(NA, nrow=nrow(cor_matrix), ncol=ncol(cor_matrix))
for (i in 1:(nrow(cor_matrix)-1)) {
  for (j in (i+1):ncol(cor_matrix)) {
    p_matrix[i,j] <- r2p(cor_matrix[i,j], n)
  }
}


# FDR calculation
# Only the upper triangular part of the matrix should be calculated, and the diagonal (autocorrelation) should be ignored.
p_matrix_adj <- p.adjust(p_matrix[upper.tri(p_matrix)], method = "fdr")


# Return adjusted p-values to the matrix
adj_matrix <- matrix(NA, nrow=nrow(cor_matrix), ncol=ncol(cor_matrix))
adj_matrix[upper.tri(adj_matrix)] <- p_matrix_adj


# Selection of gene pairs with statistical significance
threshold <- 0.05
sig_cor <- which(adj_matrix < threshold, arr.ind = TRUE)

# Result extraction
results <- data.frame(
  Gene1 = rownames(cor_matrix)[sig_cor[,1]],
  Gene2 = rownames(cor_matrix)[sig_cor[,2]],
  Correlation = cor_matrix[sig_cor],
  AdjustedPvalue = adj_matrix[sig_cor]
)

#=======================
# Save results
#=======================
# Save all results
#write.table(results, "/path/to/your/output/significant_correlations_FDR0.05.txt", quote=F, sep="\t", row.names=F)

# Save results, excluding pairs where both genes are associated with IR events
results_fil <- results[!(grepl("\\.IR", results$Gene1) & grepl("\\.IR", results$Gene2)), ]
write.table(results_fil, output_file, quote=FALSE, sep="\t", row.names=FALSE)

