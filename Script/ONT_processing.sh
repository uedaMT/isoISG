#!/bin/bash
#SBATCH --job-name=FLAIR

#==================================================
# Author: Mahoko T. Ueda
# Date: Mar. 27, 2024
# Description: This script runs FLAIR, a bioinformatics tool for analyzing nanopore sequencing data. 
# FLAIR must be installed prior to running this script. The job can be submitted to a SLURM cluster for execution.
# 
# Note: Prior to running this script, ensure that FLAIR is installed and that the necessary input data is available 
# in the specified directories. 
#=================================================


genome=""
fastq=""
junc=""
quality=1
out="flair-nanopore"
dir="out"
cpu=4
min=""
gtf=""


# Function to print help
print_help() {
  echo "Usage: flair_isoISG.sh -f <fastq1[, fastq2, ...]> [-t <cpu>] [-g <genome>] [-h]"
  echo "Options:"
  echo "  -f, --fastq              Input FASTQ file(s) (comma/space separated)"
  echo "  -G, --gtf                isoISG GTF file or your own gtf"
  echo "  -j, --junc               Junction files made by STAR (.out.tab)"  
  echo "  -c, --cage               TSS data (obtained from https://reftss.riken.jp/datafiles)"  
  echo "  -g, --genome Reference   genome file"
  echo "  -q, --quality            MAPQ cutoff (default: 1)"
  echo "  -t, --cpu                CPU number (default: 4)"
  echo "  -o, --out                OUTPUT file (default: flair-nanopore)"
  echo "  -d, --dir                OUTPUT directory (default: ./out)"  
  echo "  -h, --help               Show this help message"
  https://reftss.riken.jp/datafiles
}

# Get options
while [[ $# -gt 0 ]]; do
  case $1 in
      -f|--fastq)
      # Split comma or space separated files into an array
      IFS=', ' read -r -a files <<< "$2"
      for file in "${files[@]}"; do
        if [[ -f $file ]]; then
          fastq_files+=("$file")
        else
          echo "Invalid FASTQ file: $file"
          exit 1
        fi
      done
      shift
      shift
      ;;
	-G|--gtf)
  	  gtf="$2"
  	  shift
  	  shift
  	  ;;
	-c|--cage)
	  cag="$2"
	  shift
	  shift
	  ;;
    -g|--genome)
      genome="$2"
      shift
      shift
      ;;
    -t|--cpu)
      cpu="$2"
      shift
      shift
      ;;
    -o|--out)
      out="$2"
      shift
      shift
      ;;
    -d|--dir)
      dir="$2"
      shift
      shift
      ;;    
    -j|--junc)
      junc="$2"
      shift
      shift
      ;;
    -q|--quality)
      quality="$2"
      shift
      shift
      ;;
    -h|--help)
      print_help
      exit
      ;;
    *)
      echo "Invalid option: $1"
      exit 1
      ;;
  esac
done


# Check if required options are present
if [ ${#fastq_files[@]} -eq 0 ]; then
  echo "FASTQ files not specified"
  exit 1
fi

#if [ -z "$junc" ]; then
#  echo "Junction file not specified"
#  exit 1
#fi

if [ -z "$genome" ]; then
  echo "Genome file not specified"
  exit 1
fi


# Print options
echo "Reference genome: $genome"
echo "GTF file: $gtf"
echo "Junction file: $junc"
echo "TSS data: $cag"
echo "FASTQ files: ${fastq_files[@]}"
echo "CPU: $cpu"
echo "Quality: $quality"
echo "Output: $out"
echo "Out directoty: $dir"


#////////// Mapping  /////////////

mkdir -p ${dir}
python flair.py align -g "${genome}" \
       -r "${fastq_files[@]}" \
       -o "${dir}/${out}" \
       --quality "${quality}" \
       -t "${cpu}" -p


#////////// Junction Correction /////////////

mkdir -p ${dir}/correct
python flair.py correct -q $dir/${out}.bed \
       -g ${GENOME} \
       -c psl/hg38-${ASM}_chrom_size.txt \
       -o ${dir}/correct/${out} \
       -f ${gtf} \
       -j ${jnc} \
       -t ${cpu}


#////////// Collapse  /////////////

mkdir -p ${dir}/collapse
python flair.py collapse -g ${GENOME} \
       -r ${FASQ} \
       -q ${dir}/correct/${out}_corrected.psl \
       -p ${cag} \
       -o ${dir}/collapse/${out} \
       -f ${gtf} \
       -t ${cpu} --generate_map  --keep_intermediate

