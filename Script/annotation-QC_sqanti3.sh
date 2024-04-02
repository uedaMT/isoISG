#!/bin/bash
#SBATCH --job-name SQANTI3
#SBATCH -o %J.sqanti3.out
#SBATCH -e %J.sqanti3.err  

######################################################################
#  performing quality control using SQANTI3 
######################################################################
#
# This script performs  quality control using SQANTI3 
# Please ensure that you have the necessary dependencies installed before running this script.
# Author: Mahoko T. Ueda
# Date: March 27, 2024
#
#
###############################
# Set your parameters
###############################

cpu=10
NAM=name_of_your_study
GTF=/path/to/the/directory/gencode_referebce_annotation.gtf
GENOME=/path/to/the/directory/gencode_reference_genome.fa
GTF=/path/to/the/directory/cupcake_gff/flnc.clst.hq.collapsed.gff
#JNK=/path/to/the/directory/junction_file/STAR_SJ.out.tab
JNK=$(find /path/to/the/directory/junction_file/ -name "*_SJ.out.tab" | paste -sd "," -) #If there are multiple junction files
CAGE=/path/to/the/directory/cage_file
PLM=/path/to/the/directory/polyA_motir_list  #obtained from https://github.com/ConesaLab/SQANTI3/tree/master/data/polyA_motifs
PLP=/path/to/the/directory/polyA_peak  #obtained from https://github.com/ConesaLab/SQANTI3/tree/master/data/polyA_motifs

#===============================================================================
# activate SQANTI environment....
#export PATH=$PATH:/path/to/directory/SQANTI3/cDNA_Cupcake/sequence
#export PYTHONPATH=$PYTHONPATH:/path/to/directory/SQANTI3/cDNA_Cupcake/sequence
#===============================================================================

KEY=sample_name
mkdir -p ${NAM}

#////////////// Annotation  /////////////////
python sqanti3_qc.py ${GTF} \
       ${GTF} ${GENOME} \
       --cage_peak ${CAGE} \
       --polyA_motif_list ${PLM} \
       --polyA_peak ${PLP} \
       -o ${NAM}_SQANTI3 \
       --report pdf \
       --isoAnnotLite --genename \
       -c ${JNK} \
       -t ${cpu} \
       -d ${NAM}



#////////////// Filtering  /////////////////
mkdir -p ${NAM}/filter
python sqanti3_RulesFilter.py \
		${NAM}/${KEY}_classification.txt \
    	--isoforms ${NAM}/${KEY}_corrected.fasta \
    	--gtf ${NAM}/${KEY}_corrected.gtf \
    	--faa ${NAM}/${KEY}_corrected.faa \
		-d ${NAM}/filter \
		-o ${KEY} \
		-j /path/to/filter_custom.json

