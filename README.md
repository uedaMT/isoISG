# isoISG
Welcome to our repository, which serves as a source of annotation data, scripts, and supplementary data for our paper titled "Functional and Dynamic Profiling of Transcript Isoforms Reveals Essential Roles of Alternative Splicing in Interferon Response."

In this study, we generated a novel isoform annotation for IFN-α-stimulated and unstimulated B-cells called isoISG (isoforms of Interferon-Stimulated Genes). This was accomplished using the advanced PacBio sequencing platforms (Sequel II/IIe), enabling us to explore the complex landscape of alternative splicing events triggered by the interferon response.

## Download

To explore the isoISG annotation further, we have made available the following files:<br><br>

**The isoISG annotation files generated in LCL**
                                  
- GTF file: [isoISG.gtf.gz](https://drive.google.com/u/4/uc?id=1Yaw3TFNB3AT9HVHPWEOqVgpPsb0mm8D7&export=download)
- GTF file with strict condition [isoISG-strict.gtf.gz](https://drive.google.com/u/4/uc?id=1zkaXl88swa0I5o6oPvfn679ii0chgUfe&export=download)<br>
  This file contains isoforms that are stringently filtered to include only those with both CAGE peaks and polyA motifs.
- SQANTI classification file: [isoISG_sqanti-classification.txt.gz](https://drive.google.com/u/4/uc?id=1nIWPZVXKruxbNjZOjhn7BYiJDq7oIDFV&export=download)<br><br><br>

**Isoform sequences**
- Nucleotide sequence: [isoISG_nt.fa.gz](https://drive.google.com/u/4/uc?id=1447GPoYqbjyqlhskpcAh22UCywfieuuH&export=download)
- Amino acid sequence: [isoISG_aa.fa.gz](https://drive.google.com/u/4/uc?id=1w0BJhcenNjnMJXkOtqg_9hporbnJUJRW&export=download)<br><br><br>

**Data for DEG/DET analysis using isoISG**
- Gene and isoform table: [isoISG_txt2gene.txt.gz](https://drive.google.com/u/4/uc?id=1rnLK59YDGbGuUn4pvCc05RhjQ9wzjlhS&export=download)
- Isoforms with IR in the CDS of reference gene: [isoISG_txt2gene.txt.gz](https://drive.google.com/u/4/uc?id=1blInW5zaI_qFeUWz8es8qGqJZOoHKLtQ&export=download)
<br><br><br>

**Supplementary data for our paper**
- Data.gz
- Table.gz<br><br><br>


### Scripts for analysis
The scripts used for analysis in our study can be found in this repository. They are instrumental for replicating our analysis or adapting it to your own datasets.
- [Access the Scripts](https://github.com/uedaMT/isoISG/tree/main/Script)<br><br><br>


## Raw data availability
These raw datasets are available in the DDBJ Sequence Read Archive (DRA), and can be accessed using the provided accession numbers.
The following sequencing datasets were generated and analyzed in this study:
- **PacBio Sequel II/IIe Iso-Seq data**: Accession number: [DRA016394](https://identifiers.org/insdc.sra:DRA016394)
- **ONT RNA-seq**: Accession number: [DRA016393](https://identifiers.org/insdc.sra:DRA016393)
- **Short-read RNA-seq data**: Accession number: [DRA016395](https://identifiers.org/insdc.sra:DRA016395)<br><br><br>


## Visualization tool:
A UCSC browser is available: https://genome.ucsc.edu/s/UEDA/isoISG<br><br>

If you're interested in our research, please check out our full paper for an overview of our main findings. Thank you for visiting our repository!
