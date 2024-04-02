#!/bin/bash
#SBATCH --job-name=minimap2_featureCounts
#SBATCH -o %J.out
#SBATCH -e %J.err  

##################################################################################
# Mapping and Quantifying PacBio Iso-Seq Reads Using Minimap2 and FeatureCounts
##################################################################################
#
# This script performs the following operations for each sample:
# 1. Aligns PacBio Iso-Seq reads to a reference genome using Minimap2.
# 2. Sorts the alignments using Samtools.
# 3. Quantifies exon-level expression using FeatureCounts.
#
# Ensure that the necessary software (minimap2, samtools and FeatureCounts) are installed.
# The reference genome (for Minimap2) and the GTF annotation file (for FeatureCounts)
# must be specified correctly.
#
# Usage:
# Adjust the 'CPU', 'REF', 'GTF', and paths to input/output directories as needed.
#
# Author: Mahoko T. Ueda
# Date: March 27
#
##################################################################################

CPU=10 # Update CPU number for your environment
REF="/path/to/your/reference_genome.fasta"
GTF="/path/to/your/isoISG.gtf"
SAMPLE_LIST="/path/to/your/sample_list.txt"
OUTPUT_DIR="/path/to/your/output"


# Ensure output directories exist
mkdir -p ${OUTPUT_DIR}

while read LINE; do
    CCS_READS=/path/to/your/CCS_read/${LINE}.hifi_reads.fasta.gz  # CCS reads in compressed FASTA format
    BAM="${OUTPUT_DIR}/${LINE}_ccs-mapped-hg38.bam"
    
    # Align CCS reads with Minimap2 and sort with Samtools
    minimap2 -ax map-pb -t ${CPU} --secondary=no ${REF} ${CCS_READS} \
    | samtools view -Su - \
    | samtools sort -@ ${CPU} -O BAM -o ${BAM}
    
    # Quantify with FeatureCounts
    featureCounts -T ${CPU} \
                  -t exon \
                  -a ${GTF} \
                  -g transcript_id \
                  -o ${OUTPUT_DIR}/${LINE}_counts_exon.txt ${BAM}

done < "${SAMPLE_LIST}"
