#############################################
#    Perform IsoformSwitchAnalyzeR analysis
#############################################
#
# This script performs DIU and functional consequence analyses using IsoformSwitchAnalyzeR.
# Please ensure that you have the necessary dependencies installed before running this script.
# Author: Mahoko T. Ueda
# Date: March 27, 2024
#
# Usage: Rscript IsoformSwitchAnalyzeR.R [options]
#        ex) Rscript IsoformSwitchAnalyzeR.R --help
#
#############################################

args <- commandArgs(trailingOnly = TRUE)

# Define default values
dir <- getwd()
input <- NULL
output <-NULL

# Parse the options
for (i in seq_along(args)) {
  if (args[i] == "--dir") {
    dir <- args[i + 1]
  } else if (args[i] == "--input") {
    input <- args[i + 1]
  } else if (args[i] == "--output") {
    output <- args[i + 1]
  } else if (args[i] == "--gtf") {
    gtf_file <- args[i + 1]
  } else if (args[i] == "--fasta") {
    fasta_file <- args[i + 1]
  } else if (args[i] == "--list") {
    file_list <- args[i + 1]
  } else if (args[i] == "--help") {
    cat("Usage: Rscript IsoformSwitchAnalyzeR.R [options]\n")
    cat("\nOptions:\n")
    cat("  --dir DIR: Working directory\n")
    cat("  --input DIRNAME: Directory where your quant.sf files (output of Salmon) are located.\n")
    cat("  --output FILENAME: Output file name  (e.g., [output]_AS_table_withName.txt)\n")
    cat("  --gtf FILENAME: GTF file name\n")
    cat("  --fasta FILENAME: Fasta file name\n")
    cat("  --list FILENAME: Table containing sample name and group\n")
    q(save = "no")
  }
}

# Use the options
parentDir = input
setwd(dir)

#=======================================
# File explanations
#=======================================
# 1) Table containing sample name and group (--list)
#  	 Design file (sampleID condition) 
#	 sampleID	condition
#	 SRR7733637	unstimulated
#	 SRR7733639	unstimulated
#	 SRR7733640	unstimulated
#	 SRR7733606	24h
#	 SRR7733607	24h
#	 SRR7733608	24h
#
# 2) Use the line of code below to convert the 'condition' column into a factor variable with custom levels.
#	 In this case, the custom levels are 'unstimulated', '24h', and '72h'.
#    myDesign$condition<- factor(myDesign$condition, levels = c("unstimulated","24h","72h"))
#
# 3) Prepare isoform quantification matrix by salmon.
#
# 4) Gene name for isoISG (isoISG_txt2gene.txt), or you can create your own dataset.
#	 TXNAME 	GENEID
#	 PB.2.1	ENSG00000276256
#	 PB.2.2	ENSG00000276256
#	 PB.2.3	ENSG00000276256
#	 PB.3.1	novelGene_7
#	 PB.3.4	novelGene_9	
#
# 5) GTF file (isoISG_annot.gtf) or your own gtf file.
#
# 6) fasta file (isoISG_nt.fa) or your own fasta file.
#
# 7) Othre requirements
# 	 Prediction of protein domains ===> Pfam can be run either locally or via their webserver.  
# 	 Prediction of coding Potential ===> CPAT can be run either locally or via their webserver.
# 	 Prediction of Signal Peptides ===> SignalP can be run either locally or via their webserver (V5 is supported)
# 	 Prediction of Intrinsically Disordered Regions ===> IUPred2A can run either locally or via their webserver
#========================================


######################################
#         Script Starts Here
######################################
#install.packages("IsoformSwitchAnalyzeR")
library(IsoformSwitchAnalyzeR)


#--------------------------
# Read salmon data
#--------------------------
# Read the Salmon quantification data by importing the isoform expression.
# Provide the 'parentDir' variable as the path to the directory where your quant.sf files are located.

salmonQuant <- importIsoformExpression(parentDir,     
	showProgress = FALSE,
    ignoreAfterBar = FALSE,
    ignoreAfterPeriod = FALSE)


#--------------------------
# Imporeting data
#--------------------------
# Import the necessary data for IsoformSwitchAnalyzeR analysis.
# Provide the required parameters such as 'isoformCountMatrix', 'isoformRepExpression', 'designMatrix',
# 'isoformExonAnnoation', and 'isoformNtFasta'.
# Modify the paths to the GTF file and the FASTA file according to your specific file locations.


# Set the paths to the GTF and fasta files
gtf_path <- file.path(dir, gtf_file)
fasta_path <- file.path(dir, fasta_file)



SwitchList <- importRdata(
    isoformCountMatrix   = salmonQuant$counts,
    isoformRepExpression = salmonQuant$abundance,
    designMatrix         = myDesign,
    isoformExonAnnoation = gtf_path,
    isoformNtFasta       = fasta_path,
    showProgress = FALSE,
    addAnnotatedORFs = TRUE,
    ignoreAfterBar = FALSE,
    ignoreAfterPeriod = FALSE,
	ignoreAfterSpace = FALSE
)

#--------------------------
# Filtering data
#--------------------------

SwitchList <- preFilter(
    switchAnalyzeRlist = SwitchList,
    geneExpressionCutoff = 1,
    isoformExpressionCutoff = 0,
    removeSingleIsoformGenes = TRUE
)


#--------------------------
#  Perform DEXseq analysis
#--------------------------
# Perform DEXseq analysis using the IsoformSwitchAnalyzeR package.
# Provide the required parameters such as 'switchAnalyzeRlist', 'alpha', 'dIFcutoff', and others.
# Modify the parameters based on your specific analysis requirements.

SwitchList <- isoformSwitchTestDEXSeq(
		#--- Core arguments ---
        switchAnalyzeRlist = SwitchList,
        alpha = 0.05,
        dIFcutoff = 0.05,
        #--- Advanced arguments ---
        correctForConfoundingFactors=TRUE,  #e.g. batch effects, Default = TRUE
        overwriteIFvalues=TRUE,	            #Correct batch effect corrected IF and dIF, Default = TRUE
        reduceToSwitchingGenes = TRUE,      #Works on dIF values corrected for confounding effects if overwriteIFvalues=TRUE. Default is TRUE.
        reduceFurtherToGenesWithConsequencePotential = TRUE,
        onlySigIsoforms = FALSE,    #====> Both isoforms the pairs considered if reduceFurtherToGenesWithConsequenc should be significantly differential used. Default = FALSE
        keepIsoformInAllConditions = TRUE,  #Default = TRUE
        showProgress = TRUE,   #Default = FALSE
        quiet = FALSE
)


# Save the resulting SwitchList object as an RDS file.
# Modify the 'file' parameter to specify the desired output file name and path.
# saveRDS(SwitchList, file = "dexseq_result_dIF0.05_p0.05.rds")


#------------------------------------------
#  Add info for the significant sequences
#------------------------------------------

SwitchList <- analyzeCPAT(
    switchAnalyzeRlist   = SwitchList,
    pathToCPATresultFile = "/path/to/your/cpat_prediction",
    codingCutoff         = 0.725,
    removeNoncodinORFs   = TRUE
)

SwitchList <- analyzeSignalP(
    switchAnalyzeRlist   = SwitchList,
    pathToSignalPresultFile = "/path/to/your/signalP_prediction"
)

SwitchList <- analyzeIUPred2A(
    switchAnalyzeRlist   = SwitchList,
	pathToIUPred2AresultFile="/path/to/your/iupred2a_prediction",
	showProgress = FALSE
)

SwitchList <- analyzePFAM(
    switchAnalyzeRlist   = SwitchList,
    pathToPFAMresultFile = "/path/to/your/hmm_prediction",
    showProgress=FALSE
)

#------------------------------------------
# Alternative splicing events
#------------------------------------------
# Analyze alternative splicing events using the IsoformSwitchAnalyzeR package.
# Perform alternative splicing analysis on the SwitchList object.
# Provide the required parameters such as 'SwitchList', 'alpha', and 'dIFcutoff'.

SwitchList <- analyzeAlternativeSplicing(
    SwitchList,
	alpha = 0.05,
	dIFcutoff = 0.05
)


#----------------------------------------
#  Testing switch consequence
#----------------------------------------
# Analyze switch consequences of interest using the IsoformSwitchAnalyzeR package.
# Define the consequences of interest that you want to analyze.
# Modify the 'consequencesOfInterest' vector to include the desired consequences.


consequencesOfInterest <- c('intron_retention',
'coding_potential',
'NMD_status',
'domains_identified', 
'tss', 
'tts', 
'last_exon', 
'5_utr_seq_similarity',
'5_utr_length',
'3_utr_seq_similarity',
'3_utr_seq_similarity',
'3_utr_length', 
'IDR_identified', 
'IDR_length', 
'IDR_type', 
'signal_peptide_identified')


# Analyze the switch consequences of interest on the SwitchList object.
# Provide the required parameters such as 'switchAnalyzeRlist', 'alpha', 'consequencesToAnalyze', and 'dIFcutoff'.
# Modify the parameters based on your specific analysis requirements.

SwitchList <- analyzeSwitchConsequences(
    switchAnalyzeRlist = SwitchList,
	alpha = 0.05, 
    consequencesToAnalyze = consequencesOfInterest, 
    dIFcutoff = 0.05, # very high cutoff for fast runtimes - you should use the default (0.1)
    showProgress=FALSE
)

# Save the resulting SwitchList object as an RDS file.
# Modify the 'file' parameter to specify the desired output file name and path.
saveRDS(SwitchList, file = "SwitchAnalyzeR_result_dIF0.05_p0.05.rds")

#------------------------------------------
# Adding gene names and output tables
#------------------------------------------
# Read the gene names from the "gene_name_text" file.
name <- read.table("txt2gene_isoISG.txt", header=T, sep="\t")
colnames (name) <- c('isoform_id', 'gene_name')
names (SwitchList)

DIU <- SwitchList$isoformSwitchAnalysis
DIU <- merge (DIU, name, by="isoform_id")
dd <- subset (DIU, DIU$padj < 0.05)

DOM <- SwitchList$domainAnalysis
DOM <- merge (DOM, name, by="isoform_id")

CON <- SwitchList$switchConsequence #already with gene_name 
cc <- subset (CON, CON$switchConsequence != 'NA')

AS <- SwitchList$AlternativeSplicingAnalysis
AS <- merge (AS, name, by="isoform_id")

#head (CON$switchConsequence)
#write.table (DIU, "AS-event_table_Switch.txt", quot = F, row.names=F,sep="\t")
#write.table (DOM, "AS-event_table_domain.txt", quot = F, row.names=F,sep="\t")
#write.table (CON, "AS-event_table_consequence.txt", quot = F, row.names=F,sep="\t")
#Output slaced TPM matrix
#write.table (SwitchList$isoformRepExpression, "scaled_TPM.txt", quot = F, row.names=F, sep="\t")


write.table(dd, paste0(output, "_DIU_table_pval0.05.txt"), quot = FALSE, row.names = FALSE, sep = "\t")
write.table(cc, paste0(output, "_Consequence_TRUE.txt"), quot = FALSE, row.names = FALSE, sep = "\t")
write.table(AS, paste0(output, "_AS_table_withName.txt"), quot = FALSE, row.names = FALSE, sep = "\t")



######################################
#         Script Ends Here
######################################
