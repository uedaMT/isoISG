#############################################
#    Perform DEG analysis
#############################################
#
# This script performs DEG analysis using DESeq2.
# Please ensure that you have the necessary dependencies installed before running this script.
# Author: Mahoko T. Ueda
# Date: May 30, 2023
#
# Usage: Rscript DESeq2.R [options]
#        ex) Rscript DESeq2.R --help
#
#############################################


#-----------------------------------------------------------
#  Files and info you need to run this script
#-----------------------------------------------------------
# 1) Prepare sample group list (tab-delniated, --list).
# 	 (path/to/your/file/sample_group.txt)
#
#    Sample	Group
#	 NA18943	control
#	 NA18944	contorl
#	 NA18943_IFNa2	IFN
#	 NA18944_IFNa2	IFN
#
# 2) Directorey where your salmon output are locate (--salmon)
#
# 3) Gene and transcript table (isoISG_txt2gene.txt) for gene level DEG analysis
#
#-----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

# Define default values
dir <- getwd()

# Parse the options
for (i in seq_along(args)) {
  if (args[i] == "--dir") {
    dir <- args[i + 1]
  } else if (args[i] == "--salmon") {
    salmon <- args[i + 1]
  } else if (args[i] == "--output") {
    output <- args[i + 1]
  } else if (args[i] == "--list") {
    file_list <- args[i + 1]
  } else if (args[i] == "--help") {
    cat("Usage: Rscript DESeq2.R [options]\n")
    cat("\nOptions:\n")
    cat("  --dir DIR: Working directory\n")
    cat("  --salmon FILENAME: Directory where your quant.sf files (output of Salmon) are located. \n")
    cat("  --output FILENAME: Output file name\n")
    cat("  --list FILENAME: Table containing sample name and group\n")
    q(save = "no")
  }
}


######################################
#         Script Starts Here
######################################

library("tximport")
library("DESeq2")

# Use the options
setwd(dir)

samples <- read.table(file_list, header = T, sep="\t")
files <- file.path(salmon, samples$Sample, "quant.sf")

# Read tx2gene mapping
tx2gene <- read.table("isoISG_txt2gene.txt", header=T, sep="\t")


#-- gene level comparison
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

#-- transcript level comparison
#txi <- tximport(files, type = "salmon", txOut = TRUE)


col_num <- ncol(txi$counts)
countData <- round (as.matrix(txi$counts))
colData <- read.table(file_list,row.names=1, header=T)
colnames(countData) <- row.names(colData)


if (!all(rownames(colData) %in% colnames(countData))) {
  stop("Mismatch between sample names in the group list and Salmon output files")
}
countData <- countData[, rownames(colData)]
if (!all(rownames(colData) == colnames(countData))) {
  stop("Mismatch between sample names in the group list and count data")
}

# Run DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData,colData, design= ~Group)
dds <- DESeq(dds)
keep <- rowSums (counts (dds)) >=col_num
dds <- dds [keep,]
# Get normalized counts
# c <- counts(dds, normalized=TRUE)
# write.table(c, paste0(args$output_name, "_normalized_counts.txt"), quote=F, sep="\t")
 
res <- results(dds)
res <- res[order(res$padj), ]
res1 <- as.data.frame(res)
res2 <- subset (res1, res1$padj < 0.05)

write.table (res2, paste0(output,"_FDR0.05.txt"), quote=F, sep="\t")
write.table (res1,  paste0(output,"_all.txt"), quote=F, sep="\t")
#write.table (res2, "DESeq2_DE-trx_FDR0.05.txt"), quote=F, sep="\t")
#write.table (res1,  "DESeq2_DE-trx_all.txt"), quote=F, sep="\t")

