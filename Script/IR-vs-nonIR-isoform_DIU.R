#################################################################
#  Perform Differential IR- vs nonIR-Isoform usage analysis 
#################################################################
#
# This script performs DIU analyses (IR vs non-IR).
# Please ensure that you have the necessary dependencies installed before running this script.
# Author: Mahoko T. Ueda
# Date: May 30, 2023
#
#=======================================
# Usage:
#=======================================
# 1) Set the working directory to the location where your data files are stored.
#    To do this, replace the path within the setwd() function below with the desired directory path.
	 setwd("")

# 2) Prepare group file (group_IFN.txt) 
	# make text file of the condition group list (group_IFN.txt) like following.
	# gr
	#    sample_id condition
	#1  SRR7733637   control
	#2  SRR7733639   control
	#3  SRR7733640   control
	#4  SRR7733606       24h
	#5  SRR7733607       24h
	#6  SRR7733608       24h
	#7  SRR7733609       24h
	#8  SRR7733621       72h
	#9  SRR7733622       72h
	#10 SRR7733623       72h
	#11 SRR7733624       72h

# 3) Specify the type and make groups for comparison.
	gr <- read.table ("group_IFN.txt", header=T)	
	#-- Make group Control vs 24h
	grp <- gr[c(1:7),]
	typ <- "24h"

	#-- Make group Control vs 72h
	#grp <- gr[c(1:3,8:11),]
	#typ <- "72h"
	
#========================================



######################################
#         Script Starts Here
######################################

library(data.table)
library(GenomicFeatures)
library(tximport)
library(DRIMSeq)
library(stageR)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggbeeswarm)


# Read raw count data for IR and nonIR isoforms
c <- read.table ("raw_cnt_isoform_IR-vs-nonIR.txt", header=T)
cts <- c[,c(3:13)]
rownames (cts) <- c$ID
cou <- cts[,colnames(cts) %in% grp$sample_id]
txgn <- c[,c(1,2)]
colnames (txgn) <- c("GENEID", "TXNAME")

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
#min_samps_gene_expr: Minimal number of samples where genes should be expressed. 
#min_samps_feature_expr: Minimal number of samples where features should be expressed.
#min_samps_feature_prop: Minimal number of samples where features should be expressed.
#min_gene_expr: Minimal gene expression.
#min_feature_expr: Minimal feature expression.
#min_feature_prop: Minimal proportion for feature expression. This value should be between 0 and 1.
#run_gene_twice: Whether to re-run the gene-level filter after the feature-level filters.
#minor_allele_freq: Minimal number of samples where each of the genotypes has to be present.
#		  
#n = nrow(gr)
#n.small = min(table(gr$condition))
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
write.table (res1, "drimseq_result_IR_FDR0.05_cont-", typ, "_minGeneCount10_miniSmpGene3_miniSmplPro0.01.txt", quote=F, sep = "\t", row.names=F)

# All result
#write.table (res1, "drimseq_result_IR_all_cont-", typ, "_minGeneCount10_miniSmpGene3_miniSmplPro0.01.txt", quote=F, sep = "\t", row.names=F)


