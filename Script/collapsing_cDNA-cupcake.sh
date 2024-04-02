#!/bin/bash
#SBATCH --job-name cupcake
##SBATCH --job-name minimap2
###SBATCH -o %J.out
####SBATCH -e %J.err  

######################################################################
#  mapping Ribo-seq data on the human genome
######################################################################
#
# This script performs cupcake cllapsing
# Please ensure that you have the necessary dependencies installed before running this script.
# Author: Mahoko T. Ueda
# Date: Mar 27, 2024
#
#
###############################
# Set your parameters
###############################

CPU=10
# High-quality fasta of isoseq3 cluster output
FAST=/path/to/the/directory/[NAME].flnc.clst.hq.fasta  #change the file name
BAM=/path/to/the/directory/[NAME].flnc.clst.hq.bam  #change the file name


#================================== Activate cupcake ================================
#export PATH=$PATH:/path/to/the/directory/cDNA_Cupcake/sequence
#export PYTHONPATH=$PYTHONPATH:/path/to/the/directory/cDNA_Cupcake/sequence
#====================================================================================

mkdir -p cupcake
collapse_isoforms_by_sam.py --input ${FAST} \
	       --max_5_dif 100 \
	       --max_3_diff 50 \
	       --bam ${BAM} -o cupcake/flnc.clst.hq
