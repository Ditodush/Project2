
## Authors


Maha Abdelrhman, Dimitri Dushuashvili

# SARS-CoV-2 genome assembly from Illumina & ONT reads (project 2)

## Overview


In the following project our goal is to construct a genome sequence using a reference-based approach. For this purpose we will go thrue the basic tools and commands. It will give as appartunity to process Illumina and Nanopore sequencing data.

By the end of the project we will generate a SARS-CoV-2 consensus genome file from FASTQ raw sequencing data. 

The following tools will be used. We will discuss them in detail.

### QC

- FastQC (Illumina)
- fastp (Illumina)
- NanoPlot (Nanopore)
- Filtlong (Nanopore)
- BAMclipper (Illumina & Nanopore)

### Mapping

- minimap2 (Illumina & Nanopore)
- SAMtools (Illumina & Nanopore)
- IGV (Illumina & Nanopore)

### Variant calling & Consensus

- Medaka (Nanopore)
- freebayes (Illumina)
- BCFtools (Illumina & Nanopore)
 
### Lineage annotation

- Pangolin (Illumina & Nanopore)


We decided to use mamba as an updated version of conda. It will help to speed up the process. 


```
# at the first step we config some channels
conda config --add channels default
conda config --add channels bioconda
conda config --add channels conda-forge

# now we create a new environment and only install mamba in it and a specific python version 
conda create -n workshop mamba python=3.9
# then we activate the environment
conda activate workshop
# ... and use mamba to install all tools
mamba install fastqc fastp nanoplot pycoqc filtlong minimap2 samtools bcftools igv pangolin president snpeff bamclipper freebayes

```

## Input data

```
- Illumina

- Nanopore

```


### With a following command we set now variables.

```
ILLUMINA_SAMPLE1='illumina.R1.fastq.gz'
ILLUMINA_SAMPLE2='illumina.R2.fastq.gz'
NANOPORE_SAMPLE='nanopore.fastq.gz'
```

## FASTQ quality control

```
# activate the conda environment
conda activate workshop

fastqc -t 4 $ILLUMINA_SAMPLE1 $ILLUMINA_SAMPLE2
```

After the above command, we can check the HTML output file.

Next our goal is to quality-trim the Illumina reads to get rid of low-quality base calls especially at the end of reads.

### Quality trimming

```
fastp -i $ILLUMINA_SAMPLE1 -I $ILLUMINA_SAMPLE2 -o clean_reads.R1.fastq.gz -O clean_reads.R2.fastq.gz --thread 4 --qualified_quality_phred 20 

fastqc -t 2 clean_reads.R{1,2}.fastq.gz
```


### Nanopore

- NanoPlot

#### Quality assessment

```
# activate the conda environment
conda activate workshop

# run NanoPlot on your FASTQ file
NanoPlot -t 4 --fastq $NANOPORE_SAMPLE -o nanoplot/raw 
    
# run NanoPlot on your FASTQ file with some more parameters
NanoPlot -t 4 --fastq $NANOPORE_SAMPLE --title "Raw reads" \
    --color darkslategrey --N50 --loglength -o nanoplot/raw 
```

###Length filtering

- filtlong

```
# Here we need to write wy we change a length .... 
filtlong --min_length 300 --max_length 600 $NANOPORE_SAMPLE | gzip - > clean_reads_nanopore.fastq.gz

NanoPlot -t 4 --fastq clean_reads_nanopore.fastq.gz --title "Filtered reads" \
    --color darkslategrey --N50 --loglength -o nanoplot/clean 
```
