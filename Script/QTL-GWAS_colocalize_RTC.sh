#!/bin/bash
#SBATCH --job-name=QTLtools

######################################################################
#  Perform RTC analysis to check colocalzation between GWAS and sQTL
######################################################################
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
DIR2=/path/to/your/vcf_file/expression_file

VCF=VCF_file
BED=Expression_file
QTL=QLT_outout
GWAS=QTLtools_ext_all_associations_v1.0.2.1_hg38_sort.txt
HOT=hotspots.bed #(lifted-to-hg38: hotspots_b37_hg19.bed provided by QTLtools )
KEY=GWAS #study name

###############################


mkdir -p out

for i in $(seq 1 22); do
#for i in $(seq 4 7); do
    QTLtools rtc --vcf ${DIR1}/chr${i}_${VCF}.gz \
	 --bed ${DIR2}/chr${i}_${BED}.gz \
	 --hotspot ${HOT} \
	 --gwas-cis ${GWAS} \
	  ${QTL} \
	  --normal \
	  --out out/chr${i}_rtc_${KEY}.txt
	
	  # Merge results
      awk '{ if ($20 >= 0.9){print $0}}' out/chr${i}_rtc_${KEY}.txt >out/chr${i}_rtc_${KEY}_0.9.txt
done

# Extract the RTC score >0.9
cat out/chr*_rtc_GWAS_${KEY}_0.9.txt >  all_significant_rtc_GWAS_${KEY}_0.9.txt

