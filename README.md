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
|---|---|---|---|
| ($1) OUTDIR | Directory you want read counts to be outputted to | | required |
| ($2) BAM_DIR | Directory containing the .bam files (minus /align/) | | required |
| ($3) FRAC_COUNTS | Whether you want reads mapping to multiple locations to be counted as a fractional amount. (If not these counts will be discarded) | 'TRUE' or 'FALSE' | required |
| ($4) STRANDEDNESS | The strandedness of the RNAseq fastq data | unstranded=0, forward=1, reverse=2 | required |
| ($5) GTFFILE | GTF file for reference genome | | required |
| ($6) THREADS | Number of threads to run script on | | required |

**Example**
```bash
read_counts.sh "GSE282218/counts" "out_dir" FALSE 2 "ref/ncbi_dataset/data/GCF_000001405.40/genomic.gtf" 16
```
Output will be a read_counts.txt and a read_counts.txt.summary file.

For the final script in the biological processing, you will need to switch environments
```bash
conda deactivate
conda activate environment_r
```

6. Run `create_colmap_gse282218.r`
This is the one script that gets custom made for each new dataset. This one is hard coded specifically to process the GSE282218 data. \
If making one of your own create_colmap scripts, it must have a few key elements for the data frame formatting output.
- The number of rows must equal the number of columns in `read_counts.txt`
- The first column must be named "accession" and the names in the column must match perfectly the headers of `read_counts.txt`
- all coefficients one wishes to examine/compare must be under a column named coef_# (exp. coef1)
- The reference category one wishes to compare their coefficients to must be under the column ref_coef# (exp. ref_coef1)
- their must be a 'replicate' column numbering the replicates for each coefficient subtype.
- any other columns for ones own meta data management are fine to be present so long as they do not have one of the aforementioned column names.
- the table must be tab separated

**Example**
```bash
./create_colmap_gse282218.r --count_dir "GSE282218/counts" --meta_dir "GSE282218/meta" --out_dir "GSE282218/meta"
```
Output will be a `colData.tsv` file.

8. Run `differential_expression_analysis.R`
This script uses DeSEQ2 to run the differential expression analysis of the read count profile.
(DEGs = Differentially Expressed Genes)

| Argument | Description | Options | Default |
|---|---|---|---|
| --count_file | read count profile (`read_counts.txt`) | | required |
| --out_dir | output directory for DEG files | | "./out" |
| --nskip | Number of commented lines to skip before getting to the column header line. (Should run on default) | | 1 |
| --alpha | Significance cutoff for DEGs. P-values and False Discovery Rate | | 0.05 |
| --col_map | File containing meta data about column headers + samples | | required |

```bash
./differential_expression_analysis.R --count_file "GSE282218/counts/read_counts.txt" --out_dir "GSE282218/counts/de_out" --col_map "GSE282218/meta/colData.tsv"
```
Output files will include: \
- A gene of interest file (used for ORA and GeneSICA) --> `ORA_coef1_vs_control_genes_of_interest.txt`
- An ORA gene reference set (used for ORA) --> `ORA_reference_genes.txt`
- A GSEA rank file --> `GSEA_coef1_vs_Control_genes_of_interest.rnk`
- A text file containing extra haplotypes which were removed --> `removed_alt_haps_assmbl.txt`
- An extra data folder containing the full results (including insignificant) and meta data about parameters used. -->  `extra_data/full_results_coef1_vs_Control.txt` and `extra_data/result_params_coef1_vs_Control.txt` 

### Silver Simulated Dataset Setup


### GeneSICA Method
This section will describe how to run the GeneSICA method on the prepared biological and Silver datasets.

First, activate the python environment.
```bash
conda deactivate
conda activate environment_python
```

A) To run the GeneSICA + SIC Ordering method, run `boostrap_SIC_parallel.py`.

| Argument | Description | Options | Default |
|---|---|---|---|
| --goi_file | Text file containing genes of interest that are differentially expressed | | required |
| --gmt | GMT file, gene set database | | required |
| -m | Size of resample DEG list, L' | An integer | 55 |
| -n | Number of bootstrap iterations | An integer | 1000 |
| --minsize | Minimum size of gene set to be considered | An integer | 5 |
| --maxsize | Maximum size of gene set to be considered | An integer | 500 |
| --out_dir | File containing meta data about column headers + samples |  | "out/pvalues" |

**Example**
```bash
cd genesica
./bootstrap_SIC_parallel.py \
   --goi_file "../biological_data/GSE282218/counts/de_out/ORA_coef1_vs_control_genes_of_interest.txt" \
   --gmt "../ref/gmt/c5.go.v2026.1.Hs.symbols.gmt" \
   -m 55 \
   -n 100000 \
   --n_jobs 40 \
   --minsize 5 \
   --maxsize 500 \
   --out_dir "results"
```

B) To run the GeneSICA + GPD method, run `boostrap_SIC_feed_gpd.py`.

| Argument | Description | Options | Default |
|---|---|---|---|
| --goi_file | Text file containing genes of interest that are differentially expressed | | required |
| --gmt | GMT file, gene set database | | required |
| -m | Size of resample DEG list, L' | Integer > 0 | 55 |
| -n | Number of bootstrap iterations | Integer > 0 | 1000 |
| --minsize | Minimum size of gene set to be considered | An integer | 5 |
| --maxsize | Maximum size of gene set to be considered | An integer | 500 |
| --out_dir | File containing meta data about column headers + samples |  | "out/pvalues" |
| --max_tail | Number of permutated SIC values to start GPD goodness of fit test | An Integer | 500 |
| --transformation_exp | Exponent used to transform the GPD. To get a heavier tailed GPD and more reliable results, ideally you want a shape parameter < 0.  | Int > 0 | 1 |

**Example**
```bash
cd genesica
./bootstrap_SIC_parallel.py \
   --goi_file "../biological_data/GSE282218/counts/de_out/ORA_coef1_vs_control_genes_of_interest.txt" \
   --gmt "../ref/gmt/c5.go.v2026.1.Hs.symbols.gmt" \
   -m 55 \
   -n 100000 \
   --n_jobs 40 \
   --minsize 5 \
   --maxsize 500 \
   --out_dir "results" \
   --transformation_exp 1
```

### Alternative Gene Set Analysis Methods
This section explains how to run alternative gene set enrichment methods: ORA and GSEA.

First, activate the R environment
```bash
conda deactivate
conda activate environment_r
```

1. Run `webgestalt.R`
| Argument | Description | Options | Default |
|---|---|---|---|
| --method | The gene set enrichment method |'ORA' or 'GSEA'| 'ORA' |
| --goi_file | Text file containing genes of interest that are differentially expressed | | required |
| --goi_type | The ID format that the proteins in the goi_file are in | | ensembl_gene_id |
| --db_file | A gmt file containing gene set database info |  | required |
| --enrich_db_type | The ID format that the proteins in the db_file are in  | | ensembl_gene_id |
| --reference_file | A text file containing the reference set of genes. Only required for ORA | | |
| --reference_type | The ID format that the proteins in the reference_file are in | | ensembl_gene_id |
| --out_dir | Output directory for the results |  | "webgestalt_out" |
| --project_name | Project name used in directory and file naming |  | "" |
| --minNum | Minimum size of gene set to be considered  | Int > 0 | 5 |
| --maxNum| Maximum size of gene set to be considered  | Int > 0 | 500 |
| --perNum | Number of permutations for null distribution. GSEA only  | Int > 0 | 1000 |
| --threads | Number of threads to run script on  | Int > 0 | 1 |



## Authors
Johanna Enright, Katie Ovens, Farhad Maleki.

## Citations
1. Signal, B., & Kahlke, T. (2022). how_are_we_stranded_here: quick determination of RNA-Seq strandedness. BMC bioinformatics, 23(1), 49. https://doi.org/10.1186/s12859-022-04572-7
