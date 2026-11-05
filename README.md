# GeneSICA
Scripts used to run the GeneSICA algorithm as well as data processing scripts, alternative gene set analysis scripts, and plotting scripts.

## Overview
GeneSICA aims to identify biologically relevant enrichment in smaller gene sets with less genes overlap. The method combines an information-theoretic score with weighted bootstrapping to assess statistical significance.
This is because large gene sets containing genes shared across many biological processes may provide limited insight into the specific processes underlying phenotypic differences in case-control studies. Other methods such as Over-representation analysis make simplifying assumptions, including that all genes are equally important when performing gene set analysis, which we hope to address with this method.

## Installation
There are three different conda environments used in the various processing stages of this project.
1. `environment_bio_processing.yml` is used when downloading RNA-seq data, cleaning fasta files, running genome alignments, and performing read counts.
2. `environment_r.yml` is used when running any file ending in .r
3. `environment_python.yml` is used when running any file ending in .py 


```bash
git clone https://github.com/JohannaEnright/GeneSICA.git
cd geneSICA
conda env create -f environment_bio_processing.yml
conda env create -f environment_r.yml
conda env create -f environment_python.yml
```

## Usage

### Data Download and Processing
All scripts to download and process RNA-seq data can be found in the `biological_data` folder. 

activate the environment_bio_processing environment
```bash
conda activate environment_bio_processing
```

1. run the `download_geo_experiment.sh` script. \
This script only processes paired-end reads. Output files can be quite large so it is recommended to download them on the scratch if available.

| Argument | Description | Default |
|---|---|---|
| ($1) GSE | GEO accession number | required |
| ($2) SRP | SRP | required |
| ($3) OUTDIR | Output directory | required |
| ($4) THREADS | Number of Threads | required |

**Example**
```bash
bash download_geo_experiment.sh "GSE282218" "SRP546017" "out_dir" 4
```
The output is paired-end read zipped fastq files in `out_dir/fastq/` and separate meta data files at `out_dir/meta/srr_list.txt` and `out_dir/meta/${SRP}metadata.tsv`

2. Run `clean_fastqs.sh`
This script ensures all fastqs in the file are intact and performs quality filtering using default parameters.

| Argument | Description | Default |
|---|---|---|
| ($1) INPUT_DIR | Directory where the fastq files are located (minus the `fastq` dir) | required |

**Example**
```bash
bash clean_fastqs.sh "out_dir/fastq"
```
The output will be directed to `out_dir/cleaned`. JSON and HTML files summarizing cleaning for each fastq will also be outputted.

3. Run `star_build_index.sh`
This script builds the index necessary for read alignment.

| Argument | Description | Default |
|---|---|---|
| ($1) FASTA | Fasta files from the reference genome of the organism of choice | required |
| ($2) GTF | GTF files for reference genome of choice. Headers much match those in the fasta | required |
| ($3) SJDB_OVERHANG | max_read_length - 1 (100 works if unknown) | required |
| ($4) SCRATCH | Output directory for index | required |
| ($5) THREADS | Number of Threads to run script on | required |

**Example**
```bash
bash star_build_index.sh "ref/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna" "ref/ncbi_dataset/data/GCF_000001405.40/genomic.gtf" 92 out_dir 16
```
The output will be directed to `out_dir/star_index/${SJDB_OVERHANG}/`

4. run `star_align_human.sh`
This script maps the cleaned fastq reads to their genomic locations using the index built in step 3.

| Argument | Description | Default |
|---|---|---|
| ($1) OUTDIR | Directory that cleaned fastq and meta files are in | required |
| ($2) STAR_INDEX | File path that the genomic index file is located in | required |
| ($3) THREADS | Number of threads to run script on | required |

**Example**
```bash
star_align_human.sh "out_dir" "out_dir/star_index/92" 16
```
Outputs .bam and .csi files to `out_dir/align`

4.5 *Optional: Run `check_strand.sh`*
This script uses the python package how_are_we_stranded_here^1 to determine the strandedness of RNAseq fastq data.
If the strandedness of your data is unknown, follow instructions for download and use here: [how_are_we_stranded_here](https://github.com/signalbash/how_are_we_stranded_here).

This information is recommended to have for step 5.

5. Run `read_counts.sh`

| Argument | Description | Options | Default |
|---|---|---|
| ($1) OUTDIR | Directory you want read counts to be outputted to | | required |
| ($2) BAM_DIR | Directory containing the .bam files (minus /align/) | | required |
| ($3) FRAC_COUNTS | Whether you want reads mapping to multiple locations to be counted as a fractional amount. (If not these counts will be discarded) | 'TRUE' or 'FALSE' | required |
| ($4) STRANDEDNESS | The strandedness of the RNAseq fastq data | unstranded=0, forward=1, reverse=2 | required |
| ($5) THREADS | Number of threads to run script on | | required |

**Example**
```bash
read_counts.sh "GSE282218/counts" "out_dir" FALSE 2 16
```
Output will be a read_counts.txt and a read_counts.txt.summary file.


## Authors
Johanna Enright, Katie Ovens, Farhad Maleki.

## Citations
1. Signal, B., & Kahlke, T. (2022). how_are_we_stranded_here: quick determination of RNA-Seq strandedness. BMC bioinformatics, 23(1), 49. https://doi.org/10.1186/s12859-022-04572-7
