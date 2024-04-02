#!/bin/bash
#SBATCH --job-name PacBio
#SBATCH -o %J.out
#SBATCH -e %J.err


######################################################################
# Processing PacBio Iso-Seq Data
######################################################################
#
# This script processes PacBio Iso-Seq data for subsequent analysis.
# Steps include demultiplexing with lima and refining with isoseq3 refine.
#
# Adjust the 'LINE' variable to specify the sample of interest.
#
# Ensure that the necessary software (lima and isoseq3) are installed and
# accessible within the conda environment 'isoseq3'.
#
# Author: Mahoko T. Ueda
# Date: March 27
#
######################################################################

CPU=10 # Update CPU number for your environment

# Path to the sample list file
FILENAME="/path/to/sample_list.txt"

# Adjusted base directory for processing to point to PacBio CCS reads
DIR_BASE="/path/to/your/PacBio_ccs_read"

while read LINE; do
    DIR="${DIR_BASE}/${LINE}"

    # Create a directory for lima output and run lima for demultiplexing
    mkdir -p  "${DIR}/lima"
    lima --isoseq "${DIR}/TKDpj_${LINE}.hifi_reads.bam" primer.fa "${DIR}/lima/${LINE}.fl.bam" -j ${CPU}

    # Specify the full-length non-chimeric (FLNC) bam file
    FL="${DIR}/lima/${LINE}.fl.IsoSeq_5p--IsoSeq_3p.bam"

    # Create a directory for isoseq3 output and run isoseq3 refine for polishing
    mkdir -p  "${DIR}/isoseq3"
    isoseq3 refine --require-polya ${FL} primer.fa "${DIR}/isoseq3/${LINE}_refine-flnc.bam" -j ${CPU}
    
    # Add the clustering step
    mkdir -p "${DIR}/isoseq3/cluster_hq"
    isoseq3 cluster ${DIR}/isoseq3/${LINE}_refine-flnc.bam \
        ${DIR}/isoseq3/cluster_hq/${LINE}_cluster-flnc.bam \
        --verbose -j ${CPU}

done < "$FILENAME"