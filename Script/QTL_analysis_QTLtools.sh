#!/bin/bash
#SBATCH --job-name=QTLtools
#SBATCH -o stdout.%J
#SBATCH -e stderr.%J 

################################
#  Perform sQTL analysis
################################
#
# This script performs RTC analysis.
# Please ensure that you have the necessary dependencies installed before running this script.
# Author: Mahoko T. Ueda
# Date: May 30, 2023
#
###############################
# Set your parameters
###############################

DIR1=/path/to/your/vcf_file
DIR2=/path/to/your/expression_file
VCF=VCF_file
BED=Expression_file
WIN=1000000

###############################


for i in $(seq 1 22); do
#for i in $(seq 1 3); do
    mkdir -p perm
    QTLtools cis --vcf  ${DIR1}/chr${i}_${VCF}.gz \
	    --bed ${DIR2}/chr${i}_${BED}.gz \
	    --permute 1000 --window ${WIN} --out perm/chr${i}_perm1000_win${WIN}.txt --seed 123456
done


# Merge results
cat perm/chr*_perm1000_win${WIN}.txt | gzip -c > permutations_full.txt.gz

# FDR using qtltools_runFDR_cis.R (providede from QTLtools site)
qtltools_runFDR_cis.R permutations_full.txt.gz  0.05 results_FDR0.05.txt


