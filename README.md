# SLRanger

<a href="https://pypi.python.org/pypi/SLRanger" rel="pypi">![PyPI](https://img.shields.io/pypi/v/SLRanger?color=green) </a>

An integrated approach for spliced leader detection and operon prediction in eukaryotes using long RNA reads 

## Workflow

<div align="center">
  <img src="document/workflow.png" width="700" alt="Workflow">
</div>

 **(A)** A flow chart illustrates the direct RNA sequencing procedure of the species with the spliced leader and trans-splicing mechanism. SL RNA joins to pre-mRNA through trans-splicing, adding an SL sequence to the 5′ end. This forms a mature mRNA with the SL sequence, which is then sequenced using Nanopore sequencing. Long RNA reads were mapped to genome reference. 
 **(B)** A flow chart illustrates the workflow of the function of SL detection using SLRanger. After long RNA reads are mapped to reference, the 5’ unaligned end sequences will be extracted for alignment to SL sequence references and random sequences (as control). Based on the scoring system produced by SLRanger, we obtained the “SL score” on the possibility that one read has an SL sequence and the distribution of different SL types. 
 **(C)** The plot of the principle on how SLRanger predicts operon. Based on the information of whether each read has a high-confident SL sequence, we predict the operon structure through the mapping information of each read and the relative position of the gene. A high proportion of reads of the SL1 type will be regarded as operon upstream genes, while a high proportion of SL2 type reads or multiple reads supported by SL2 types will be regarded as operon downstream genes.

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
  -b BAM, --bam BAM     input the bam file (required)
  -o OUTPUT, --output OUTPUT
                        output file (default: SLRanger.txt)
  --visualization       Turn on the visualization mode
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
will be output in a webpage format. An example is provided [here](sample/SLRange_view/visualization_results.md).

####  Example
We provided test data to run as below.
```
git clone https://github.com/lrslab/SLRanger.git
cd sample/
unzip data.zip
SL_detect.py --ref SL_list_cel.fa --bam test.bam -o SLRanger.txt -t 4 --visualization
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
operon_predict.py -g cel_wormbase.gff -b test.bam -i SLRanger.txt  -o test.gff
```
