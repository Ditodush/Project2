
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






