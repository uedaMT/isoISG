#!/bin/bash
#SBATCH --job-name Salmon
#SBATCH -o %J.out
#SBATCH -e %J.err  


#############################################
# Run RNA-seq quantification using salmon
#############################################
#
# This script performs RNA-seq quantification with salmon.
# Please ensure that you have the necessary dependencies installed before running this script.
# Author: Mahoko T. Ueda
# Date: May 30, 2023
#
#############################################


#-------------------------------------------------------------
# Set the following parameters
#-------------------------------------------------------------

cpu=10
INX=/set/your/path/to/salmon_index
DIR=/set/your/path/to/fastq
LIST=/set/your/path/to/sample_list
REF=/set/your/path/to/fasta_file #(isoISG_nt.fa)
KEY=name_of_study 

#-------------------------------------------------------------

#//// Generate index ////
salmon index -t ${REF} -i ${INX}/${KEY} -p ${cpu}



#//// Quantification ////
FILENAME=${LIST}
mkdir -p salmon/${KEY}

while read LINE;do
    echo "${LINE}:"
    salmon quant -i ${INX}/${KEY} -l A \
	   -1 ${DIR}/${LINE}_1_trim.fastq.gz \
	   -2 ${DIR}/${LINE}_2_trim.fastq.gz \
	   -p ${cpu} --validateMappings \
	   --seqBias --gcBias \
	   --writeUnmappedNames \
	   -o salmon/${KEY}/${LINE}
done < $FILENAME



#/////////////////  This is the end of the program... ////////////////////

