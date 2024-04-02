######################################################################
# Running ribomap for RNA and ribosome profiling data integration
######################################################################
#
# This script generates a CDS range file from a given FASTA file, which is
# used in integrating RNA-seq and ribosome profiling data using ribomap (ribomap_sc.sh).
# It assumes the entire transcript is coding and calculates the length accordingly.
# Please ensure that you have the necessary dependencies installed before running this script,
# including Biopython.
#
# Author: Mahoko T. Ueda
# Date: March 27, 2024
#
#
###############################
# Set your parameters
###############################


from Bio import SeqIO

def create_cds_range_file(fasta_file, output_file):
    """
    Generates a CDS range file from a FASTA file.
    
    Parameters:
    fasta_file (str): Path to the input FASTA file.
    output_file (str): Path to the output file for storing CDS ranges.
    """
    with open(output_file, 'w') as out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            transcript_id = record.id
            cds_length = len(record.seq)
            # start position is 0, end position is length-1 (0-based indexing)
            out.write(f"{transcript_id} 1 {cds_length}\n")


# Input and output file paths
fasta_file = "/path/to/your/fasta/your_data.fa"  # Update this path to your fasta file
output_file = "/path/to/your/output/cds_ranges.txt"          # Update this path for the output file

create_cds_range_file(fasta_file, output_file)
