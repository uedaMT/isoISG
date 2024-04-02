#################################################################
#  Perform Differential IR- vs nonIR-Isoform usage analysis 
#################################################################
#
# This script performs DIU analyses (IR vs non-IR).
# Please ensure that you have the necessary dependencies installed before running this script.
# Author: Mahoko T. Ueda
# Date: May 30, 2023
#
# Usage: Rscript IR-vs-nonIR-isoform_DIU.R [options]
#        ex) Rscript IR-vs-nonIR-isoform_DIU.R --help
#
#################################################################

args <- commandArgs(trailingOnly = TRUE)

# Define default values
dir <- getwd()
input <- "raw_count.txt"
output <- "drimseq_result.txt"

# Parse the options
for (i in seq_along(args)) {
  if (args[i] == "--dir") {
    dir <- args[i + 1]
  } else if (args[i] == "--input") {
    input <- args[i + 1]
  } else if (args[i] == "--output") {
    output <- args[i + 1]
  } else if (args[i] == "--list") {
    file_list <- args[i + 1]
  } else if (args[i] == "--help") {
    cat("Usage: Rscript IR-vs-nonIR-isoform_DIU.R [options]\n")
    cat("\nOptions:\n")
    cat("  --dir DIR: Working directory\n")
    cat("  --input FILENAME: Input file name\n")
    cat("  --output FILENAME: Output file name\n")
    cat("  --list FILENAME: Table containing sample name and group\n")
    q(save = "no")
  }
}

# Use the options
setwd(dir)

#=======================================
# File explanation
#=======================================
# 1) Set the working directory (--dir) to the location where your data files are stored.
# 2) Specify the type or names (--typ) for comparison.	 
# 3) Prepare group file (--list) 
#
#	   sample_id condition
#	 SRR7733637   control
#	 SRR7733639   control
#	 SRR7733640   control
#	 SRR7733606       24h
#	 SRR7733607       24h
#	 SRR7733608       24h
#	 SRR7733609       24h
#
# 4) Prepare raw read matrix for IR and non-IR isoforms (--input) 
#    #Raw read counts for IR- and nonIR-isoforms
#		  Gene	ID	SRR7733637	SRR7733639	SRR7733640	SRR7733606	SRR7733607	SRR7733608	SRR7733609
#	 	  AAAS	AAAS.IR	91.239	181.155	26.264	85.132	153.435	58.614	49.533
#		  AAAS	AAAS.nonIR	307.761	259.846	164.736	274.707	223.565	281.386	448.462
#		  AASDHPPT	AASDHPPT.IR	53.259	47.984	100.932	102.425	60.517	20.799	98.141
#		  AASDHPPT	AASDHPPT.nonIR	773.197	688.017	458.328	1051.577	899.485	629.204	1109.874
#		  ABCC5	ABCC5.IR	0	118.143	26.093	0	46.599	68.943	181.637
#		  ABCC5	ABCC5.nonIR	334	184.857	119.907	198.001	157.4	116.056	209.363	
#========================================



######################################
#         Script Starts Here
######################################
library(data.table)
library(DRIMSeq)
library(dplyr)


# Read the comparison file specified by --list
if (!is.null(file_list)) {
  grp <- read.table(file_list, header = TRUE)  # Update the function and arguments based on your file format
}


# Read raw count data for IR and nonIR isoforms
#cts <- read.table ("raw_cnt_isoform_IR-vs-nonIR.txt", header=T)
# Read the input file specified by --input
if (!is.null(input)) {
  cts <- read.table(input, header = TRUE)  # Update the function and arguments based on your file format
  txgn <- cts[,c(1,2)]
  rownames (cts) <- cts$ID
  cou <- cts[,colnames(cts) %in% grp$sample_id]
  colnames (txgn) <- c("GENEID", "TXNAME")
}


dim(cou)
range(colSums(cou)/1e6)

txgn.sub = txgn[match(rownames(cou),txgn$TXNAME),]
counts = data.frame(gene_id = txgn.sub$GENEID, feature_id = txgn.sub$TXNAME, cou)

d <- dmDSdata(counts = counts, samples = grp)
methods(class=class(d))
#counts(d[1,])[,1:4]


#/////////////////////////////////////////////////////////////////////////////////////////
#  First filter the object, before running procedures to estimate model parameters
#/////////////////////////////////////////////////////////////////////////////////////////
#
# (1) it has a count of at least 10 in at least n.small samples, 
# (2) it has a relative abundance proportion of at least 0.1 in at least n.small samples,
# (3) the total count of the corresponding gene is at least 10 in all n samples.
#
# min_samps_gene_expr: Minimal number of samples where genes should be expressed. 
# min_samps_feature_expr: Minimal number of samples where features should be expressed.
# min_samps_feature_prop: Minimal number of samples where features should be expressed.
# min_gene_expr: Minimal gene expression.
# min_feature_expr: Minimal feature expression.
# min_feature_prop: Minimal proportion for feature expression. This value should be between 0 and 1.
# run_gene_twice: Whether to re-run the gene-level filter after the feature-level filters.
# minor_allele_freq: Minimal number of samples where each of the genotypes has to be present.
#		  
# n = nrow(gr)
# n.small = min(table(gr$condition))
#-------------------------------------------------------------------------------------------------


d <- dmFilter(d,
           min_samps_feature_prop=0.01, 
           min_samps_gene_expr=3, 
		   min_gene_expr=10)

design_full <- model.matrix(~ condition, data = samples(d))


#///////////////////////////////////////////////////
# Estimate the model parameters and test for DTU.
#///////////////////////////////////////////////////

set.seed(123)
system.time({
	d <- dmPrecision(d, design = design_full, BPPARAM=BiocParallel::MulticoreParam(20))
	d <- dmFit(d, design = design_full, verbose = 1, BPPARAM=BiocParallel::MulticoreParam(20))
})

contrast = colnames(design_full)[2]
set.seed(123)
system.time({
	d <- dmTest(d, coef = contrast, verbose = 1, BPPARAM=BiocParallel::MulticoreParam(20))
})

res <- results(d,level="feature")
table(res$adj_pvalue < 0.05)
res1 <- subset (res, res$adj_pvalue < 0.05)

# only significant result (FDR < 0.05)
write.table (res1, out, quote=F, sep = "\t", row.names=F)



