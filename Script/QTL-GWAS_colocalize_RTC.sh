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

#study name
KEY=GWAS 

#Directories
DIR1=/path/to/your/vcf_file
DIR2=/path/to/your/expression_file
OUT=/path/to/your/output

# Files
VCF=/path/to/your/VCF_file
BED=/path/to/your/Expression_file
QTL=/path/to/your/QLT_outout
GWAS=/path/to/your/GWAS_fata.txt
HOT=/path/to/your/hotspots_file

###############################


mkdir -p ${OUT} 

for i in $(seq 1 22); do
    QTLtools rtc --vcf ${DIR1}/chr${i}_${VCF}.gz \
	 --bed ${DIR2}/chr${i}_${BED}.gz \
	 --hotspot ${HOT} \
	 --gwas-cis ${GWAS} \
	  ${QTL} \
	  --normal \
	  --out ${OUT}/chr${i}_rtc_${KEY}.txt

	
	  # Merge results
      awk '{ if ($20 >= 0.9){print $0}}' out/chr${i}_rtc_${KEY}.txt >out/chr${i}_rtc_${KEY}_0.9.txt
done

# Extract the RTC score >0.9
cat out/chr*_rtc_GWAS_${KEY}_0.9.txt >  all_significant_rtc_GWAS_${KEY}_0.9.txt

