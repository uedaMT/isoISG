#!/bin/bash
#SBATCH --job-name=salmon
##SBATCH -o stdout.%J
##SBATCH -e stderr.%J 

######################################################################
# Quantifying Nanopore RNA-seq Data Using Salmon
######################################################################
#
# This script aligns Nanopore RNA-seq reads to a transcriptome reference
# using minimap2, followed by expression quantification with Salmon.
#
# Ensure that the necessary software (minimap2, samtools and salmon) are installed
#
######################################################################


CPU=10 # Update CPU number for your environment

# Path to the reference transcriptome
REF="/path/to/your/reference.fa"

# Path to the sample list
SAMPLE_LIST="/path/to/your/sample_list.txt"

# Directory for output
OUTPUT_DIR="/path/to/output"

while read CND; do
    echo "Processing ${CND}"
 
	mkdir -p ${OUTPUT_DIR}/${QS}/${CND}    
	FASQ=${DIR1}/${CND}_qs10_500bp_3101267.fastq.gz
    
	# Align Nanopore Reads with minimap2 to the transcriptome
	minimap2 -t ${CPU} -ax map-ont --secondary=no ${REF} ${FASQ} | \
	    samtools view -@ ${CPU} -bu - | \
	    samtools sort -n -@ ${CPU} -o ${OUTPUT_DIR}/${QS}/${CND}/${CND}_alignment_sort.bam -
    
	
	   
    # Quantify with Salmon using the alignment
    salmon quant --ont -t ${REF} \
           -a ${OUTPUT_DIR}/${CND}/${CND}_alignment_sort.bam \
           --noErrorModel \
           -o ${OUTPUT_DIR}/${CND} -p ${CPU} -l U

done < "${SAMPLE_LIST}"
