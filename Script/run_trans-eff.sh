#!/bin/bash
#SBATCH --job-name tran-eff

#################################################################
#  Perform translation efficiency prediction 
#################################################################
#
# This script performs translation efficiency prediction using ribosomeProfilingQC
# Please ensure that you have the necessary dependencies installed before running this script.
# This script requires [translat-efficiency.R]
# Author: Mahoko T. Ueda
# Date: May 30, 2023
#
#################################################################


for chr in $(seq 1 22); do
#for chr in $(seq 18 19); do
	echo "calculation starts for chr${chr}"
	Rscript --slave --vanilla translat-efficiency.R ${chr}
done
