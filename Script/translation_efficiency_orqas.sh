#!/bin/bash
#SBATCH --job-name ORQAS
#SBATCH -o %J.out
#SBATCH -e %J.err

######################################################################
# ORQAS Toolkit for RNA and Ribosome Profiling Data Integration
######################################################################
#
# This script utilizes the ORQAS toolkit for integrating RNA-seq and ribosome profiling data.
# Please ensure that you have the necessary dependencies installed before running this script.
# Author: Mahoko T. Ueda
# Date: March 27, 2024
#
#
###############################
# Set your parameters
###############################

# Directories and file paths
prefix= /path/to/the/directory/orqas_output
fasta_file=/path/to/your_data.fa


FILENAME=path/to/sample_list.txt

# Convert TXT to CDS
python ORQAStoolkit.py TXTtoCDS \
       -i ${fasta_file}  \
       -o ${prefix}

# Process each line in FILENAME
while read LINE; do
   # Aggregate CDS based on abundance
   python ORQAStoolkit.py aggregateCDS \
	  -i ${prefix}/${LINE}/${LINE}.Abundance.txt \
	  -c ${prefix}/${LINE}/${LINE}.ENSTtoCDS.txt \
	  -a ${fasta_file} \
	  -o ${prefix}/${LINE}/${LINE}.AggregateCDS.txt

   # Validate ORF
   python ORQAStoolkit.py validateORF \
	  -i /path/to/the/directory/orqas_output/${LINE}/ribomap_out.base \
	  -c ${prefix}/${LINE}/${LINE}.ENSTtoCDS.txt \
	  -o ${prefix}/${LINE}/${LINE}.validateORF.txt

   # Calculate OPM
   python ORQAStoolkit.py OPMcalculator \
	  -i ${prefix}/${LINE}  \
	  -o ${prefix}/${LINE}/${LINE}
done < $FILENAME
