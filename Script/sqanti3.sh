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
# Date: May 30, 2023
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
JNK=/path/to/the/directory/junction_file/STAR_SJ.out.tab
CAGE=refTSS_v3.3_human_coordinate.hg38.bed #obtained from https://reftss.riken.jp/datafiles
PLA=polyA_list  #obtained from https://github.com/ConesaLab/SQANTI3/tree/master/data/polyA_motifs

#===============================================================================
# activate SQANTI environment....
#export PATH=$PATH:/path/to/directory/SQANTI3/cDNA_Cupcake/sequence
#export PYTHONPATH=$PYTHONPATH://path/to/directory/SQANTI3/cDNA_Cupcake/sequence
#===============================================================================

KEY=gencode-all
mkdir -p ${NAM}

#////////////// Annotation  /////////////////
python sqanti3_qc.py ${GTF} \
       ${GTF} ${GENOME} \
       --cage_peak ${CAGE} \
       --polyA_motif_list data/polyA.list \
       -o ${NAM}_SQANTI3 \
       --report pdf \
       --isoAnnotLite --genename \
       -c ${JNK} \
       -t ${cpu} \
       -d ${NAM}



#////////////// Filtering  /////////////////
python sqanti3_RulesFilter.py \
       -a 0.6 -m 100 -c 3 \
       --report pdf \
       ${NAM}/${NAM}_classification.txt \
       ${NAM}/${NAM}_corrected.fasta \
       ${NAM}/${NAM}_corrected.gtf





