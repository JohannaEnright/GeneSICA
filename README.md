# GeneSICA
Scripts used to run the GeneSICA algorithm as well as data processing scripts, alternative gene set analysis scripts, and plotting scripts.


## Overview
GeneSICA aims to identify biologically relevant enrichment in smaller gene sets with less genes overlap. The method combines an information-theoretic score with weighted bootstrapping to assess statistical significance.
This is because large gene sets containing genes shared across many biological processes may provide limited insight into the specific processes underlying phenotypic differences in case-control studies. Other methods such as Over-representation analysis make simplifying assumptions, including that all genes are equally important when performing gene set analysis, which we hope to address with this method.

## Installation
```bash
git clone https://github.com/JohannaEnright/GeneSICA.git
cd geneSICA
conda env create -f environment_bio_processing.yml
conda env create -f environment_r.yml
conda env create -f environment_python.yml
conda activate environment_bio_processing
```

## Usage
```bash
Rscript scripts/run_sica.R --input data/expression.csv --components 20 --out results/
```

| Argument | Description | Default |
|---|---|---|
| `--input` | Genes × samples matrix (CSV) | required |
| `--components` | Number of components | 10 |
| `--out` | Output directory | `results/` |

## Input / Output
Describe file formats and what each output file contains.

## Example
Run on the demo data in `example/` to reproduce the figure above.

## Repository structure
```
scripts/   analysis scripts
data/      example input
results/   outputs
```

## Citation
If you use this, please cite: Author et al., Journal, Year. DOI

## License
MIT. See [LICENSE](LICENSE).

