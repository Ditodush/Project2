
# Project2

SARS-CoV-2 genome assembly from Illumina & ONT reads

- Here we describe basic tools and commands to process either Illumina or Nanopore sequencing data using SARS-CoV-2 amplification data as an example. The goal is to construct a genome sequence using a reference-based approach. Althought the methodological steps are quite similar, both Illumina and Nanopore need different tools and parameters.

```
# config some channels, this might be already done.
# basically helps to not explicitly type the channels
conda config --add channels default
conda config --add channels bioconda
conda config --add channels conda-forge
```
