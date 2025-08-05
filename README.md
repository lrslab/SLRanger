# SLRanger

<a href="https://pypi.python.org/pypi/SLRanger" rel="pypi">![PyPI](https://img.shields.io/pypi/v/SLRanger?color=green) </a>

An integrated approach for spliced leader detection and operon prediction in eukaryotes using long RNA reads 

## Workflow

<div align="center">
  <img src="document/workflow.png" width="700" alt="Workflow">
</div>

Workflow of SLRanger. (A) The direct RNA sequencing workflow for species exhibiting trans-splicing. In such species, spliced leader (SL) sequences are located at the 5′ ends of RNA reads due to the trans-splicing mechanism. Mature mRNAs containing SL sequences are sequenced using Oxford Nanopore direct RNA sequencing, and the resulting long reads are subsequently aligned to the reference genome. (B) The SL detection module of SLRanger. After aligning long RNA reads to the reference genome, the unaligned 5′-end fragments are extracted and aligned against a reference set of SL (including SL1 sequence and all SL2 variant sequences) and random sequences (as controls). Using SLRanger’s scoring scheme, each candidate read is assigned an “SLRanger score.” Reads will be assigned to the corresponding SL-type according to the highest “SLRanger score”. Reads that can’t be classified as SL1 will be considered as the SL2 type. If SL1 and SL2 can’t be distinguished, the read is labeled SL_unknown. If variants of SL2 can’t be distinguished, the read is labeled SL2_unknown. A dynamic cutoff is then applied to identify “high-confidence SL reads”. (C) The principle underlying operon prediction by SLRanger. Based on the presence of high-confidence SL sequences in each read and their genomic mapping positions relative to gene annotations, operon structures are inferred. Genes with a high proportion of SL1-type reads are predicted to be upstream operon genes, whereas genes with a high proportion of SL2-type reads or supported by multiple SL2-type reads are predicted to be downstream operon genes.

## Installation
 The pipeline is invoked using a CLI written in **Python(3.9-3.11)** and requires a Unix-based operating system. For conda method, we provided installation from pypi and git
###  i. Conda method
1. Prepare a new conda env
```
conda create -n SLRanger_env python=3.9
conda activate SLRanger_env
conda install -c bioconda bedtools minimap2 samtools
```
2a. Install from **PyPI**  
```
# install from pypi
pip install SLRanger
```
2b. Install from **Github**
```
#  install from git
git clone https://github.com/lrslab/SLRanger.git
cd SLRanger/
python setup.py install
```
###  ii. Docker method
```
docker pull zhihaguo/slranger_env
```
##  Manual 
SLRanger encompasses two primary functions, spliced leader (SL) detection and operon prediction, used to determine whether long RNA reads carry SL sequences and predict the operon structure based on the SL information.
### 1. Preprocessing
#### Reference and annotation selection 
The long RNA reads will be mapped to the genome reference. The genome reference (**fasta/fa/fna** file) and annotation file (**GFF** file) should be determined before running SLRanger.
These can be downloaded from [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/) or assembled independently.

In our sample folder,we provided _C. elegans_ annotation file.
#### Long reads alignment
Additionally, we require users to provide their own alignment file (**BAM** file). For long reads, minimap2 is the recommended software. 
In the sample folder, we have provided **test.bam**, which was generated using the following command.
```
minimap2 -ax splice -uf -t 80 -k14 --MD --secondary=no $reference $basecall_file > tmp.sam
samtools view -hbS tmp.sam | samtools sort -@ 32 -F 260 -o test.bam
samtools index test.bam
```
### 2. Spliced Leader detection
`SL_detect.py` is designed to detect spliced leaders. 
#### Command options
Available options can be viewed by running `SL_detect.py -h` in the command line.
```
SL_detect.py -h
usage: SL_detect.py [-h] -r REF -b BAM [-o OUTPUT] [--visualization] [-t CPU]

help to know spliced leader and distinguish SL1 and SL2

options:
  -h, --help            show this help message and exit
  -r REF, --ref REF     SL reference (fasta file recording SL sequence, required)
  -i BAM, --input BAM     input the bam file (required)
  -m , --mode           RNA or cDNA
  -o OUTPUT, --output OUTPUT
                        output file (default: SLRanger.txt)
  -t CPU, --cpu CPU     CPU number (default: 4)
  -c CUTOFF, --cutoff CUTOFF
                        The value used to filter high confident SL reads. 
                        The higher the value, the stricter it is. 
                        The range is between 0-10. (default: 4)                      
```
#### Output description

##### i. result table
| **col name**       | **description** |
|--------------------|-----------------|
| query_name         |Unique name of reads|
| strand             |Mapping direction of reads|
| soft_length        |Length of soft clipping at 5' end of reads (length of 5' unmapped region)|
| aligned_length     |Aligned length of reads|
| read_end           |End position of locally sequence mapped to SL sequence|
| query_length       |Length of locally sequence mapped to SL sequence|
| consensus          |Consensus sequence between the locally sequence mapped to SL sequence and SL sequence reference|
| random_sw_score    |Optimal score of 5' unmapped region mapped to random sequences obtain by Smith-Waterman algorithm (SW score)|
| random_final_score |Final score of 5' unmapped region mapped to random sequences obtained by SLRanger scoring system|
| random_SL_score    |Random final score normlized by the maximum possible score for the length of the locally mapped region sequence (SL_score)|
| sw_score           |Optimal score of 5' unmapped region mapped to SL reference sequences obtain by Smith-Waterman algorithm (SW score)|
| final_score        |Final score of 5' unmapped region mapped to SL reference sequences obtained by SLRanger scoring system|
| SL_score           |Final score normlized by the maximum possible score for the length of the locally mapped region sequence (SL_score)|
| SL_type            |Spliced Leader types; random if random_SL_score > SL_score|

##### ii. visualization result
The summary table and figures, including the Data Summary Table and the pictures including Cumulative Counts (SW), Cumulative Counts (SL), Query Length Distribution, Aligned Length Distribution, SL Type Distribution.
will be output in a webpage format. An example is provided [here](sample/SLRanger_view/visualization_results.md).

####  Example
We provided test data to run as below.
```
git clone https://github.com/lrslab/SLRanger.git
cd sample/
unzip data.zip
SL_detect.py --ref SL_list_cel.fa --input RNA_test.bam -o SLRanger.txt -t 4 
SL_detect.py --ref SL_list_cel.fa --input cDNA_test.bam -o SLRanger_cDNA.txt -t 4
```
### 3. Operon prediction
`operon_predict.py` is designed to predict operons.
#### Command options
Available options can be viewed by running `operon_predict.py -h` in the command line.
```
operon_predict.py  -h
usage: operon_predict.py [-h] -g GFF -b BAM -i INPUT [-o OUTPUT] [-d DISTANCE]
help to know spliced leader and distinguish SL1 and SL2

options:
  -h, --help            show this help message and exit
  -g GFF, --gff GFF     GFF annotation file (required)
  -b BAM, --bam BAM     bam file (required)
  -i INPUT, --input INPUT
                        input the SL detection file (result file from SL_detect.py, required)
  -o OUTPUT, --output OUTPUT
                        output operon detection file ( default: SLRanger.gff)
  -d DISTANCE, --distance DISTANCE
                        promoter scope (default: 5000)
  -c CUTOFF, --cutoff CUTOFF
                        The value used to filter high confident SL reads. 
                        The higher the value, the stricter it is. 
                        The range is between 0-10. (default: 4)                      
```
#### Output description
A GFF file will be returned.

####  Example
We provided test data to run as below (should be run after `SL_detect.py`).
```
cd sample/
operon_predict.py -g cel_wormbase.gff -b RNA_test.bam -i SLRanger.txt  -o test.gff
```
