#################################################################
#  Perform translation efficiency prediction
#################################################################
#
# This script performs translation efficiency prediction using ribosomeProfilingQC
# Please ensure that you have the necessary dependencies installed before running this script.
# Author: Mahoko T. Ueda
# Date: May 30, 2023
#
#################################################################

#---------------------------------------------------------------------------
#    Usage
#---------------------------------------------------------------------------
# Use [R_run_trans-eff.sh] for running thie script to specify chromosome number...
# Set the following parameters...!!!

	setwd ("/path/to/your/working_directory")
	DIR="/path/to/your/output_directory"
	DIR1="/path/to/your/gtf_directory"
	DIR2="/path/to/your/ribo-bam_directory" #ribo-seq data
	DIR3="/path/to/your/rna-bam_directory" #rna-seq data
	GTF="gtf_file name"
	BAM1="bam_file name" #mapped ribo-seq data
	BAM2="bam_file name" #mapped RNA-seq data

#----------------------------------------------------------------------------


######################################
#         Script Starts Here
######################################

suppressMessages(library(GenomicFeatures))
suppressMessages(library (ribosomeProfilingQC))
chr = commandArgs(trailingOnly=TRUE)[1]

#///// read gtf ///////
txdb <- makeTxDbFromGFF(file.path(DIR1, GTF), format="gff",organism="Homo sapiens")


#///// read ribo-seq ///////
RPF <- file.path(DIR2, BAM1)

#///// read rna-seq ///////
RNA <- ile.path(DIR3, BAM2)


#head (RPF)
#head (RNA)

#///// calculating efficiency  ///////
print("Calculating coverageDepth...")
cvg <- ribosomeProfilingQC::coverageDepth(RPF, RNA, txdb)

print("Calculating Translation Efficiency...")
TE90 <- ribosomeProfilingQC::translationalEfficiency(cvg, window = 90, normByLibSize=TRUE)

dir.create(DIR)
write.table(TE90$TE,paste0(DIR,"translation_efficiency_chr",chr,".txt"),sep="\t",row.names = TRUE,col.names = NA,quote = FALSE)
write.table(TE90$mRNA,paste0(DIR,"coverage_mRNA_chr",chr,".txt"),sep="\t",row.names = TRUE,col.names = NA,quote = FALSE)
write.table(TE90$RPFs,paste0(DIR,"coverage_RPFs_chr",chr,".txt"),sep="\t",row.names = TRUE,col.names = NA,quote = FALSE)

print("ALL done!")

#---- end of the program ----- 
