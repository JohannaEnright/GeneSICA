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
Output files will include: 
- A gene of interest file (used for ORA and GeneSICA) --> `ORA_coef1_vs_control_genes_of_interest.txt`
- An ORA gene reference set (used for ORA) --> `ORA_reference_genes.txt`
- A GSEA rank file --> `GSEA_coef1_vs_Control_genes_of_interest.rnk`
- A text file containing extra haplotypes which were removed --> `removed_alt_haps_assmbl.txt`
- An extra data folder containing the full results (including insignificant) and meta data about parameters used. -->  `extra_data/full_results_coef1_vs_Control.txt` and `extra_data/result_params_coef1_vs_Control.txt` 

### Silver Simulated Dataset Setup
This section describes how to create simulated gene expression profiles from Silver.

1. First you need a txt file of the gene sets you want to be differentially expressed, formatted like `mixedsize_mixedoverlap_halfup_1_3.txt`. This can be manually created, using `select_gene_sets_by_size.py` as guidance for gene sets to include based on gene set size and overlap score.

*Optional:* \
*To determine the sizes and overlap scores of the gene sets you wish to be differentially expressed in silver, run the script `select_gene_sets_by_size.py`.*

```bash
conda deactivate
conda activate environment_python
```

| Argument | Description | Options | Default |
|---|---|---|---|
| --gmt | "../ref/gmt/c5.go.v2026.1.Hs.symbols.gmt" |  | required |
| --size_range | The size range of the gene sets you are interested in. 2 integers separated by whitespace | Int Int | required |
| --overlap_range | The overlap range of the gene sets you are interested in. 2 floats separated by whitespace | Float Float | required |
| --overlap_score_file | Pre-calculated overlap score txt file. If not provided, scores will be automatically computed from scratch | | "" |
| --out_dir | output directory for DEG files | | "./" |

**Example**
```bash
./select_gene_sets_by_size.py \
   --gmt "../ref/gmt/c5.go.v2026.1.Hs.symbols.gmt" \
   --size_range 5 20  \
   --overlap_range 0 11  \
   --out_dir "../ref/size_overlap_filtering" \
```
This will output a `filtered_${gmt}_${size_range}size_${overlap_range}overlap.tsv` file and a pre-calculated overlap score file to be used for faster analysis next time called `overlap_scores_${gmt}.json`. Use the tsv file to determine which gene sets of specific sizes and overlaps you want to include in your gene_set_file.txt file.


2. Run `create_contrast_file.r`

```bash
conda deactivate
conda activate environment_r
```

| Argument | Description | Options | Default |
|---|---|---|---|
| --count_file | read count profile (`read_counts.txt`) | | required |
| --col_map | This is the `colData.tsv` file. It contains meta data about the `read_counts.txt` file| | required |
| --out_dir | output directory for DEG files | | "./" |
| --nskip | Number of commented lines to skip before getting to the column header line in `read_counts.txt`. (Should run on default) | | 1 |
| --col_map | This is the `colData.tsv` file. It contains meta data about the `read_counts.txt` file| | required |
| --filename_add | This is the for if you want to add any text to the filename (like GEO accession)| | "" |
| --contrast_vars | The 2 variables that are used as a case-control comparison (names must match those in `colData.tsv` column coef#.) Exp. TUMOUR WILDTYPE | 2 strings separated by whitespace | required |
| --control_case | Indicate which of the contrast_vars is a case and which is a control. 'd' for case, 'c' for control. Order must correspond with contrast_var order | Exp. 'd' 'c'  | 2 strings separated by whitespace | required |

**Example**
```bash
./create_contrast_file.r \
   --count_file "../biological_data/GSE282218/counts/read_counts.txt" \
   --col_map "../biological_data/GSE282218/meta/colData.tsv" \
   --out_dir "simulated_out" \
   --filename_add "GSE282218_sim" \
   --contrast_vars "RA" "Control" \
   --control_case "d" "c"
```

This will output a `${filename_add}_contrast.txt` file.

3. Run `normalize_count_xp_data.r`. This script normalizes the counts in a biological expression profile so that it is ready for input to silver.

| Argument | Description | Options | Default |
|---|---|---|---|
| --count_file | read count profile (`read_counts.txt`) | | required |
| --nskip | Number of commented lines to skip before getting to the column header line. (Should run on default) | | 1 |
| --filename_add | This is the for if you want to add any text to the filename (like GEO accession)| | "" |
| --contrast_map_file | contrast.txt file created in step 1 | | required |
| --out_dir | output directory | | "./" |

**Example**
```bash
./create_contrast_file.r \
   --count_file "../biological_data/GSE282218/counts/read_counts.txt" \
   --contrast_map_file "simulated_out/GSE282218_sim_contrast.txt" \
   --out_dir "simulated_out" \
   --filename_add "GSE282218_sim" \
```
Output file is `${filename_add}_profile.txt`

4. Run`create_de_fold_change_for_silver_v2.py`. This will create the file with differentially expressed genes from the specified gene sets required for input to silver.

```bash
conda deactivate
conda activate environment_python
```

| Argument | Description | Options | Default |
|---|---|---|---|
| --gene_set_file | tab separated file with gene sets wanting to be differentially enriched, manually created by you. Must follow format of `mixedsize_mixedoverlap_halfup_1_3.txt` | | required |
| --gmt_file | GMT file | | required |
| --profile_file | Normalized expression profile matrix file created in step 2. | | required |
| --out_dir | output directory | | "./" |

**Example**
```bash
./create_de_fold_change_for_silver_v2.py \
   --gene_set_file "mixedsize_mixedoverlap_halfup_1_3.txt" \
   --profile_file "GSE282218_sim_profile.txt" \
   --gmt_file "../ref/gmt/c5.go.v2026.1.Hs.symbols.gmt" \
   --out_dir "simulated_out" 
```
The output file will be `DE_gene_fold_change.txt`

5. You need to install silver from the Github repo: [silver](https://github.com/FarhadMaleki/silver).
Follow all the instructions for installation and environment creation. 
```bash
conda deactivate
conda activate silver
```

6. Run `generate_xp_dataset.py`.

| Argument | Description | Options | Default |
|---|---|---|---|
| --profile | The normalized gene expression profile `profile.txt` created in step 2 | | required |
| --contrast | The contrast file `contrast.txt` created in step 1  | | required |
| --gene_fc | The `DE_gene_fold_change.txt` file | | required |
| --lower_bound_col_name | Lower bound column name in `DE_gene_fold_change.txt` | | 'FCLower' |
| --upper_bound_col_name | Upper bound column name in `DE_gene_fold_change.txt` | | 'FCUpper' |
| --fc_id_col_name | Column name of gene ID in `DE_gene_fold_change.txt` | | 'ID' |
| --profile_id_col_name | Column name of genes in `profile.txt` | | 'ID' |
| --n_simulated_ctrls | Number of desired simulated control samples | Int > 0 | 20 |
| --n_simulated_cases | Number of desired simulated case samples | Int > 0 | 20 |
| --out_dir | output directory | | "./" |
| --replacement | Whether or not you want to sample with replacement | 'True' or 'False' | 'False' |

 
**Example**
```bash
./generate_xp_dataset.py \
   --profile "simulated_out/GSE282218_sim_profile.txt" \
   --contrast "simulated_out/GSE282218_sim_contrast.txt" \
   --gene_fc "simulated_out/DE_gene_fold_change.txt" \
   --n_simulated_ctrls 25 \
   --n_simulated_cases 25 \
   --out_dir "simulated_out" 
```
This will output two files. `simulated.profile.txt` and `simulated.contrast.txt`

7. Run `differential_expression_silver.r` This script will identify the genes that are differentially expressed in from the simulated silver expression profile.

```bash
conda deactivate
conda activate environment_r
```

| Argument | Description | Options | Default |
|---|---|---|---|
| --count_file | read simulated count profile (`simulated.profile.txt`) | | required |
| --contrast_file | The simulated contrast file `simulated.contrast.txt` | | required |
| --out_dir | output directory for DEG files | | "./de_out" |
| --nskip | Number of commented lines to skip before getting to the column header line. (Should run on default) | | 0 |
| --alpha | Significance cutoff for DEGs. P-values and False Discovery Rate | | 0.05 |

```bash
./differential_expression_silver.r \
   --count_file "simulated_out/simulated_profile.txt" \
   --contrast_file "simulated_out/simulated_contrast.txt" \
   --out_dir "simulated_out/de_out"
```
The output will be 3 main results files: `GSEA_genes_of_interest.rnk`, `ORA_genes_of_interest.txt`, and `ORA_reference_genes.txt`. As well as an `extra_data/full_results.txt` file.
From here, GeneSICA or ORA or GSEA can be run on the simulated data files.

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
| --method | The gene set enrichment method |'ORA' or 'GSEA' | 'ORA' |
| --goi_file | Text file containing genes of interest that are differentially expressed | | required |
| --goi_type | The ID format that the proteins in the goi_file are in | | ensembl_gene_id |
| --db_file | A gmt file containing gene set database info |  | required |
| --enrich_db_type | The ID format that the proteins in the db_file are in  | | ensembl_gene_id |
| --reference_file | A text file containing the reference set of genes. Only required for ORA | | |
| --reference_type | The ID format that the proteins in the reference_file are in | | ensembl_gene_id |
| --out_dir | Output directory for the results |  | "webgestalt_out" |
| --project_name | Project name used in directory and file naming |  | "" |
| --minNum | Minimum size of gene set to be considered  | Int > 0 | 5 |
| --maxNum | Maximum size of gene set to be considered  | Int > 0 | 500 |
| --perNum | Number of permutations for null distribution. GSEA only  | Int > 0 | 1000 |
| --threads | Number of threads to run script on  | Int > 0 | 1 |


## Authors
Johanna Enright, Katie Ovens, Farhad Maleki.

## Citations
1. Signal, B., & Kahlke, T. (2022). how_are_we_stranded_here: quick determination of RNA-Seq strandedness. BMC bioinformatics, 23(1), 49. https://doi.org/10.1186/s12859-022-04572-7

2. Maleki, F., Ovens, K., McQuillan, I., & Kusalik, A. J. (2021). Silver: Forging almost Gold Standard Datasets. Genes, 12(10), 1523. https://doi.org/10.3390/genes12101523

3. Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S., & Mesirov, J. P. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proceedings of the National Academy of Sciences of the United States of America, 102(43), 15545–15550. https://doi.org/10.1073/pnas.0506580102
