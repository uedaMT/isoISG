#############################################
#    Perform IsoformSwitchAnalyzeR analysis
#############################################
#
# This script performs DIU and functional consequence analyses using IsoformSwitchAnalyzeR.
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

# 2) Specify the path to your input file in the 'filename' argument of the 'read.table()' function below.
#  	 Design file (sampleID condition) 
     myDesign <- read.table("path/to/your/input_file.csv", header = TRUE, sep = "\t")


# 3) Use the line of code below to convert the 'condition' column into a factor variable with custom levels.
#    myDesign$condition<- factor(myDesign$condition, levels = c("unstimulated","24h","72h"))

# 4) Set the 'parentDir' variable to the path where your quant.sf files (output of Salmon) are located.
     parentDir = "path/to/your/salmon_output"


#//// Othre requirements  /////
# Prediction of protein domains ===> Pfam can be run either locally or via their webserver.  
# Prediction of coding Potential ===> CPAT can be run either locally or via their webserver.
# Prediction of Signal Peptides ===> SignalP can be run either locally or via their webserver (V5 is supported)
# Prediction of Intrinsically Disordered Regions ===> IUPred2A can run either locally or via their webserver
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

SwitchList <- importRdata(
    isoformCountMatrix   = salmonQuant$counts,
    isoformRepExpression = salmonQuant$abundance,
    designMatrix         = myDesign,
    isoformExonAnnoation = "/path/to/your/gtf",
    isoformNtFasta       = "path/to/your/fasta",
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

#saveRDS(SwitchList, file = "path/to/your/output/filter_geneOver1_noConsequence-noAS.rds")

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
        correctForConfoundingFactors=TRUE,  #===> (e.g. batch effects, Default = TRUE)
        overwriteIFvalues=TRUE,	            #===> (Correct batch effect corrected IF and dIF, Default = TRUE)    
        reduceToSwitchingGenes = TRUE,      #===> (Works on dIF values corrected for confounding effects if overwriteIFvalues=TRUE. Default is TRUE.)
        reduceFurtherToGenesWithConsequencePotential = TRUE,  #====> (Default = TRUE)
        onlySigIsoforms = FALSE,    #====> (Both isoforms the pairs considered if reduceFurtherToGenesWithConsequenc should be significantly differential used. Default = FALSE)
        keepIsoformInAllConditions = TRUE,  #====> (Default = TRUE)
        showProgress = TRUE,   #====> (Default = FALSE)
        quiet = FALSE
)


# Save the resulting SwitchList object as an RDS file.
# Modify the 'file' parameter to specify the desired output file name and path.
# saveRDS(SwitchList, file = "dexseq_result_dIF0.05_p0.05.rds")

#--------------------------
#  extract sequence  
#--------------------------
# Extract sequence information from the SwitchList object.
# Provide the required parameters such as 'SwitchList', 'alpha', 'dIFcutoff', and others.
# Modify the parameters based on your specific requirements.
# Set 'writeToFile' to TRUE if you want to write the extracted sequences to a file.
# Set 'pathToOutput' to specify the directory where the output files will be saved.
# Modify 'outputPrefix' to specify the desired prefix for the output files.

SwitchList <- extractSequence(
        SwitchList,
        genomeObject  = NULL,
        onlySwitchingGenes = FALSE,
        alpha = 0.05,
        dIFcutoff = 0.01,
        extractNTseq = TRUE,
        extractAAseq = TRUE,
        removeShortAAseq = TRUE,
        removeLongAAseq  = FALSE,
        removeORFwithStop= TRUE,
        addToSwitchAnalyzeRlist = TRUE,
        writeToFile = TRUE,
        pathToOutput =  "fasta",
        outputPrefix='IsoformSwitch_dIF0.05_p0.05',
)


#------------------------------------------
#  Add info for the significant sequences
#------------------------------------------

SwitchList <- analyzeCPAT(
    switchAnalyzeRlist   = SwitchList,
    pathToCPATresultFile = "path/to/your/cpat_prediction",
    codingCutoff         = 0.725,
    removeNoncodinORFs   = TRUE
)

SwitchList <- analyzeSignalP(
    switchAnalyzeRlist   = SwitchList,
    pathToSignalPresultFile = "path/to/your/signalP_prediction"
)

SwitchList <- analyzeIUPred2A(
    switchAnalyzeRlist   = SwitchList,
	pathToIUPred2AresultFile="path/to/your/iupred2a_prediction",
	showProgress = FALSE
)

SwitchList <- analyzePFAM(
    switchAnalyzeRlist   = SwitchList,
    pathToPFAMresultFile = "path/to/your/hmm_prediction",
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
# Modify the file path if needed. ('isoform_id', 'gene_name')

name <- read.table("gene_name_text")

colnames (name) <- c('isoform_id', 'gene_name')
names (SwitchList)

#> names (SwitchListAna)
# [1] "isoformFeatures"             "exons"                       "conditions"                  "designMatrix"               
# [5] "sourceId"                    "isoformCountMatrix"          "isoformRepExpression"        "runInfo"                    
# [9] "orfAnalysis"                 "isoformRepIF"                "ntSequence"                  "isoformSwitchAnalysis"      
##[13] "aaSequence"                  "domainAnalysis"              "idrAnalysis"                 "signalPeptideAnalysis"      
#[17] "AlternativeSplicingAnalysis" "switchConsequence"   
#head (SwitchList$domainAnalysis)
#head (SwitchList$switchConsequence)


DIU <- SwitchList$isoformSwitchAnalysis
DIU <- merge (DIU, name, by="isoform_id")
dd <- subset (DIU, DIU$padj < 0.05)

DOM <- SwitchList$domainAnalysis
DOM <- merge (DOM, name, by="isoform_id")

CON <- SwitchList$switchConsequence #<-- already with gene_name 
cc <- subset (CON, CON$switchConsequence != 'NA')

AS <- SwitchList$AlternativeSplicingAnalysis
AS <- merge (AS, name, by="isoform_id")

head (CON$switchConsequence)


#write.table (DIU, "AS-event_table_Switch.txt", quot = F, row.names=F,sep="\t")
#write.table (DOM, "AS-event_table_domain.txt", quot = F, row.names=F,sep="\t")
#write.table (CON, "AS-event_table_consequence.txt", quot = F, row.names=F,sep="\t")
#Output slaced TPM matrix
#write.table (SwitchList$isoformRepExpression, "scaled_TPM.txt", quot = F, row.names=F, sep="\t")


write.table (dd, "DIU_table_pval0.05.txt", quot = F, row.names=F,sep="\t")
write.table (cc, "Consequence_TRUE.txt", quot = F, row.names=F,sep="\t")
write.table (AS, "AS_table_withName.txt", quot = F, row.names=F,sep="\t")


######################################
#         Script Ends Here
######################################
