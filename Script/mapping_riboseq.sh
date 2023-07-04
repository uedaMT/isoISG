#!/bin/bash
#SBATCH --job-name STAR_mapping 
#SBATCH -o %J.out
#SBATCH -e %J.err

######################################################################
#  mapping Ribo-seq data on the human genome
######################################################################
#
# This script performs ribo-seq mapping
# Please ensure that you have the necessary dependencies installed before running this script.
# Author: Mahoko T. Ueda
# Date: May 30, 2023
#
#
###############################
# Set your parameters
###############################

cpu=10
KEY=name_of_your_study
INX="/path/to/your/index_file"
DIR1="/path/to/your/riboseq_file"
DIR2="/data03/ueda/translation/RNAseq/GSE61742/STAR"
DIR3="/data01/fastq/riboseq/STAR/GSE61742"
FILENAME="list_of_ribo-seq_samples"
rRNA="/path/to/tools/STAR/rRNA_contam" #rRNA transcript reference

######################################################################

#============================================
# filter rrna
#============================================
mkdir -p ${DIR1}

while read line;do
	STAR --runMode alignReads \
			--runThreadN ${cpu} \
			--genomeDir ${rRNA} \
			--readFilesIn ${DIR1}/${LINE}_trimmed.fq.gz \
			--readFilesCommand gunzip -c \
			--outSAMunmapped Within \
			--outFilterMultimapNmax 30 \
			--outFilterMultimapScoreRange 1 \
			--outFileNamePrefix ${DIR1}/${LINE}_rm_repbase_rrna.fastq \
			--outSAMattributes All \
			--outStd BAM_Unsorted \
			--outSAMtype BAM Unsorted \
			--outFilterType BySJout \
			--outReadsUnmapped Fastx \
			--outFilterScoreMin 10 \
			--alignEndsType EndToEnd > ${DIR1}/${RIBO}_repbase_rrna_comtam.bam

done  < $FILENAME

#============================================
# mapping to transcriptome
#============================================

mkdir -p ${DIR2}
while read line;do
    echo ${LINE}
    STAR --runMode alignReads \
	 	 --runThreadN ${cpu} \
	 	 --genomeDir ${INX} \
	 	 --readFilesIn ${DIR1}/${LINE}_rm_repbase_rrna.fastq \
	 	 --outFilterMultimapNmax 8 \
	 	 --alignSJoverhangMin 8 \
	 	 --alignSJDBoverhangMin 1 \
	 	 --sjdbScore 1 \
	 	 --outFilterMismatchNmax 4 \
	 	 --alignIntronMin 20 \
	 	 --alignIntronMax 1000000 \
	 	 --alignMatesGapMax 1000000 \
	 	 --outFileNamePrefix ${DIR2}/${LINE}_ \
	 	 --outSAMattributes All \
	 	 --outSAMtype BAM SortedByCoordinate \
	 	 --limitBAMsortRAM 32000000000 \
	 	 --outFilterType BySJout \
	 	 --outReadsUnmapped Fastx

	samtools index -@ ${cpu} ${DIR2}/${LINE}_Aligned.sortedByCoord.out.bam

done  < $FILENAME

