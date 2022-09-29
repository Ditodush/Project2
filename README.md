
## Authors


Maha Abdelrhman, Dimitri Dushuashvili

# SARS-CoV-2 genome assembly from Illumina & ONT reads (project 2)

## Overview


In the following project our goal is to construct a genome sequence using a reference-based approach. For this purpose we will go through the basic tools and commands. It will give as appartunity to process Illumina and Nanopore sequencing data.

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

### Mapping & visualization

#### Illumina & Nanopore
- minimap2
- SAMtools
- IGV



SARS-CoV-2 reference sequence will be used in FASTA format to map against. File was downloaded from the following link

https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3



There are a tools (e.g. BWA-MEM) for the short-read alignment, which may perform better then minimap2, but we decide to use only one tool for the both, Illuminda and Nanopore data.


```
# map the filtered reads to the reference genome

# Illumina
minimap2 -x sr -t 4 -a -o minimap2-illumina.sam reference.fasta clean_reads.R1.fastq.gz clean_reads.R2.fastq.gz

# Nanopore
minimap2 -x map-ont -t 4 -a -o minimap2-nanopore.sam reference.fasta clean_reads_nanopore.fastq.gz
```


Different parameter settings were used for Illumina and Nanopore data.

### Process mapping results


```
# convert mapping results (SAM format) into a binary format 
# the binary format can be faster read by a machine but not by a human
# we will use the binary format (called BAM) for visualization of the mapping

# we will use a FOR loop to process both, Illumina and Nanopore data

for SAM in minimap2-illumina.sam minimap2-nanopore.sam; do

    # get the basename of the input SAM
    BN=$(basename $SAM .sam)

    # 1) convert SAM to BAM --> check the file size after convert!
    samtools view -bS $SAM > $BN.bam

    # 2) sort the BAM
    samtools sort $BN.bam > $BN.sorted.bam

    # 3) index the BAM
    samtools index $BN.sorted.bam
done
```

### Look at the mapped reads

```
# start the Integrative Genomics Viewer (IGV)
igv &
```


Reference FASTA was loaded (that was used for the mapping) and the sorted BAM file (with the mapped reads) into IGV to look at the results.

- First, load the nCoV-2019.reference.fasta reference FASTA via “Genomes” > “Load Genome from File”
- Second, load the sorted BAM file via “File” > “Load from File”



## Primer clipping
## Need to add text !!! 
#### Illumina

```
# First, we download the primer BED scheme for Cleanplex scheme that was used
# Change to another BED file if needed!
wget --no-check-certificate https://


# It's important that the FASTA header of the reference genome 
# and the IDs in the BED file match, let's check:
head nCoV-2019.reference.fasta
head cleanplex.amplicons.bedpe

# we can see: they dont match! 
# In the reference FASTA: 'MN908947.3'
# In the BED file: 'NC_045512.2'
# So we need to replace the ID in the BED file, e.g. via
sed 's/NC_045512.2/MN908947.3/g' cleanplex.amplicons.bedpe > cleanplex-corrected.amplicons.bedpe

# check again
head nCoV-2019.reference.fasta
head cleanplex-corrected.amplicons.bedpe

bamclipper.sh -b minimap2-illumina.sorted.bam -p cleanplex-corrected.amplicons.bedpe -n 4
```

#### Nanopore

```
# First, we download the primer BED scheme for the ARTIC V1200 scheme
# Change to another BED file if needed!
wget https:/

# It's important that the FASTA header of the reference genome 
# and the IDs in the BED file match, let's check:
head nCoV-2019.reference.fasta
head nCoV-2019.bed

# now we convert this BED file into a BEDPE file needed by BAMclipper.
# The Illumina BED file we used above was already in the correct BEDPE format.
# we download a python script to do so:
wget https://
# and run it
python primerbed2bedpe.py nCoV-2019.bed --forward_identifier _LEFT --reverse_identifier _RIGHT -o nCoV-2019.bedpe

# now we can use BAMclipper - finally
bamclipper.sh -b minimap2-nanopore.sorted.bam -p nCoV-2019.bedpe -n 4
```


## Variant calling

## Need to add text !!! 
#### Illumina

```
# first re-calulate the index for the reference FASTA (sometimes issues occur bc/ of different samtool versions used)
samtools faidx reference.fasta

# now variant calling
freebayes -f nCoV-2019.reference.fasta --min-alternate-count 10 \
--min-alternate-fraction 0.1 --min-coverage 20 --pooled-continuous \
--haplotype-length -1 minimap2-illumina.sorted.primerclipped.bam > freebayes-illumina.vcf
```

#### Nanopore
```
# first, we create a new env named 'medaka' and install 'mamba' and a specific version of python needed by 'medaka'
conda create -n medaka mamba python=3.6

# we activate the env
conda activate medaka

# and install medaka
mamba install -c bioconda medaka
```

#### Call variants with Medaka


```
# first generate a file with information about potential variants
# considering the used basecalling model. You should use the matching
# model from your Guppy basecalling settings!
medaka consensus --model r941_min_hac_g507 --threads 4 --chunk_len 800 --chunk_ovlp 400 minimap2-nanopore.sorted.primerclipped.bam medaka-nanopore.consensus.hdf

# actually call the variants
medaka variant reference.fasta medaka-nanopore.consensus.hdf medaka-nanopore.vcf

# annotate VCF with read depth info etc. so we can filter it
medaka tools annotate medaka-nanopore.vcf nCoV-2019.reference.fasta minimap2-nanopore.sorted.primerclipped.bam medaka-nanopore.annotate.vcf
```

## Consensus generation

#### Illumina & Nanopore

```
# switch to the workshop env if you are not already on it
conda activate workshop

# compress the annotated VCF file (needed for the next steps)
bgzip -f medaka-nanopore.annotate.vcf
 
# index a TAB-delimited genome position file in bgz format 
# and create an index file
tabix -f -p vcf medaka-nanopore.annotate.vcf.gz

# generate the consensus
bcftools consensus -f reference.fasta medaka-nanopore.annotate.vcf.gz -o consensus-nanopore.fasta

# rename the consensus FASTA, right now the FASTA ID is still the reference
sed -i 's/MN908947.3/Consensus-Nanopore/g' consensus-nanopore.fasta
```

