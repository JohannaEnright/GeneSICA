#!/bin/bash

date

CHECK_STRANDEDNESS=check_strandedness

### Genome Assembly Reference Files ###
GTF_FILE="../ref/ensembl_genomes/Homo_sapiens.GRCh38.116.gtf"
FASTA_FILE="../ref/ensembl_genomes/Homo_sapiens.GRCh38.cdna.all.fa"    #"ref/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"
KALLISTO_IDX="../ref/kallisto_index"

### Reads to Test ###

### GSE282218 ###
#R1="/scratch/44447709/fastq/SRR31385777_1.fastq.gz"
#R2="/scratch/44447709/fastq/SRR31385777_2.fastq.gz"

#R1="/scratch/44447709/fastq/SRR31385860_1.fastq.gz"
#R2="/scratch/44447709/fastq/SRR31385860_2.fastq.gz"

### GSE300130 ###
#R1="/scratch/44629397/fastq/SRR34028021_1.fastq.gz "
#R2="/scratch/44629397/fastq/SRR34028021_2.fastq.gz"

#R1="/scratch/44629397/cleaned/SRR34028047_1_cleaned.fastq.gz"
#R2="/scratch/44629397/cleaned/SRR34028047_2_cleaned.fastq.gz"

### GSE318560 ###
R1="/scratch/44793696/fastq/SRR37123810_1.fastq.gz"
R2="/scratch/44793696/fastq/SRR37123810_2.fastq.gz"

### GSE218037 ###
#R1="/scratch/44629829/fastq/SRR22301814_1.fastq.gz"
#R2="/scratch/44629829/fastq/SRR22301814_2.fastq.gz"

#R1="/scratch/44629829/cleaned/SRR22301808_1_cleaned.fastq.gz"
#R2="/scratch/44629829/cleaned/SRR22301808_2_cleaned.fastq.gz"

### GSE265795 ###
#R1="/scratch/44629831/fastq/SRR28797710_1.fastq.gz"
#R2="/scratch/44629831/fastq/SRR28797710_2.fastq.gz"

#R1="/scratch/44629831/fastq/SRR28797720_1.fastq.gz"
#R2="/scratch/44629831/fastq/SRR28797720_2.fastq.gz"

### GSE130605
#R1="/scratch/44589797/fastq/SRR9002173_1.fastq.gz"
#R2="/scratch/44589797/fastq/SRR9002173_2.fastq.gz"

#R1="/scratch/44589797/fastq/SRR9002165_1.fastq.gz"
#R2="/scratch/44589797/fastq/SRR9002165_2.fastq.gz"


$CHECK_STRANDEDNESS \
    -g $GTF_FILE  \
    -fa $FASTA_FILE  \
    -r1 $R1  \
    -r2 $R2  \
    -k $KALLISTO_IDX \
    -p


