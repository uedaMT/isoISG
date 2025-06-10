# isoISG
Welcome to our repository, which serves as a source of annotation data, scripts, and supplementary data for our paper titled "Functional and Dynamic Profiling of Transcript Isoforms Reveals Essential Roles of Alternative Splicing in Interferon Response."

In this study, we generated a novel isoform annotation for IFN-α-stimulated and unstimulated B-cells called isoISG (isoforms of Interferon-Stimulated Genes). This was accomplished using the advanced PacBio sequencing platforms (Sequel II/IIe), enabling us to explore the complex landscape of alternative splicing events triggered by the interferon response.

## Download

To explore the isoISG annotation further, we have made available the following files:<br><br>

**The isoISG annotation files generated in LCL**
                                  
- GTF file: [isoISG.gtf.gz](https://zenodo.org/records/13282235/files/isoISG.gtf.gz?download=1)
- GTF file with strict condition [isoISG-strict.gtf.gz](https://zenodo.org/records/13282235/files/isoISG-strict.gtf.gz?download=1)<br>
  This file contains isoforms that are stringently filtered to include only those with both CAGE peaks and polyA motifs.
- SQANTI classification file: [isoISG_sqanti-classification.txt.gz](https://zenodo.org/records/13282235/files/isoISG_sqanti-classification.txt.gz?download=1)<br><br><br>

**Isoform sequences**
- Nucleotide sequence: [isoISG_nt.fa.gz](https://zenodo.org/records/13282235/files/isoISG_nt.fa.gz?download=1)
- Amino acid sequence: [isoISG_aa.fa.gz](https://zenodo.org/records/13282235/files/isoISG_aa.fa.gz?download=1)<br><br><br>

**Data for DEG/DET analysis using isoISG**
- Gene and isoform table: [isoISG_txt2gene.txt.gz](https://zenodo.org/records/13282235/files/isoISG_txt2gene.txt.gz?download=1)
- Isoforms with IR in the CDS of reference gene: [isoISG_txt2gene_IR-isoform_cds.txt.gz](https://zenodo.org/records/13282235/files/isoISG_txt2gene_IR-isoform_cds.txt.gz?download=1)
<br><br><br>


### Scripts for analysis
The scripts used for analysis in our study can be found in this repository. They are instrumental for replicating our analysis or adapting it to your own datasets.
- [Access the Scripts](https://github.com/uedaMT/isoISG/tree/main/Script)<br><br><br>


## Raw data availability
These raw datasets are available in the DDBJ Sequence Read Archive (DRA), and can be accessed using the provided accession numbers.
The following sequencing datasets were generated and analyzed in this study:
- **PacBio Sequel II/IIe Iso-Seq data**: Accession number: [DRA016394](https://humandbs.dbcls.jp/en/hum0312-v1#iso)
- **ONT RNA-seq**: Accession number: [DRA016393](https://humandbs.dbcls.jp/en/hum0312-v1#DRA016393)
- **Short-read RNA-seq data**: Accession number: [DRA016395](https://humandbs.dbcls.jp/en/hum0312-v1#DRA016395)<br><br><br>


## Visualization tool
A UCSC browser is available: https://genome.ucsc.edu/s/UEDA/isoISG<br><br>


## How to cite
If you're interested in our research, please check out our full paper for an overview of our main findings. 
- Ueda MT et al. **Functional and dynamic profiling of transcript isoforms reveals essential roles of alternative splicing in interferon response**, *Cell Genomics* (2024).　[Read the full paper]([https://doi.org/10.1016/j.xgen.2024.100345](https://www.sciencedirect.com/science/article/pii/S2666979X24002659?via%3Dihub)

