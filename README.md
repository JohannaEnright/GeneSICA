# GeneSICA
Scripts used to run the GeneSICA algorithm as well as data processing scripts, alternative gene set analysis scripts, and plotting scripts.


## Overview
GeneSICA aims to identify biologically relevant enrichment in smaller gene sets with less genes overlap. The method combines an information-theoretic score with weighted bootstrapping to assess statistical significance.
This is because large gene sets containing genes shared across many biological processes may provide limited insight into the specific processes underlying phenotypic differences in case-control studies. Other methods such as Over-representation analysis make simplifying assumptions, including that all genes are equally important when performing gene set analysis, which we hope to address with this method.

## Installation
There are three different conda environments used in the various processing stages of this project. 
i) environment_bio_processing.yml is used when downloading RNA-seq data, cleaning fasta files, running genome alignments, and performing read counts.
ii) environment_r.yml is used when running any file ending in .r
iii) environment_python.yml is used when running any file ending in .py


```bash
git clone https://github.com/JohannaEnright/GeneSICA.git
cd geneSICA
conda env create -f environment_bio_processing.yml
conda env create -f environment_r.yml
conda env create -f environment_python.yml
```

## Usage

### Data Download and Processing
1. All scripts to download and process RNA-seq data can be found in the `biological_data` folder. 

i) activate the environment_bio_processing environment
```bash
conda activate environment_bio_processing
```

ii) run the download_geo_experiment.sh script. This script only processes paired-end reads. Output files can be quite large so it is recommended to download them on the scratch.

| Argument | Description | Default |
|---|---|---|
| $1 | GEO accession number | required |
| $2 | SRP | required |
| $3 | Output directory | required |
| $4 | Number of Threads | required |

Example
```bash
bash download_geo_experiment.sh "GSE282218" "SRP546017" out_dir 4
```

The output is paired-end read zipped fastq files in out_dir/fastq/ and separate meta data files at out_dir/meta/srr_list.txt and out_dir/meta/${SRP}metadata.tsv


## Example
Run on the demo data in `example/` to reproduce the figure above.

## Authors
Johanna Enright, Katie Ovens, Farhad Maleki.


