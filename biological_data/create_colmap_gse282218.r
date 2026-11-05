##Creating the table for colData required for DESEQ2

#This script will be used to create the colData dataframe required for deSEQ2
## Creating col_map for GSE282218

## libraries ##
library(tidyverse)
library(stringr)
library(glue)



parser = ArgumentParser(description = "Creating ColData.tsv necessary for DE analysis")

# optional arguments with defaults
parser$add_argument("--out_dir",
                    default = "",  
                    type = "character",
                    help = "Output Directory for Files"
)
parser$add_argument("--meta_dir",
                    default = "",  
                    type = "character",
                    help = "directory that SRP54017.metadata.tsv is in"
)
parser$add_argument("--counts_dir",
                    default = "",  
                    type = "character",
                    help = "directory that read_counts.txt is in"
)

args = parser$parse_args()

out_dir = args$out_dir
meta_dir = args$meta_dir
count_dir = args$count_dir


############# GLOBAL VARIABLES #############
OUT_DIR = out_dir
META_TSV = file.path(meta_dir, "SRP546017.metadata.tsv")
COUNT_FILE = file.path(count_dir, "read_counts.txt")
NSKIP = 1

############################################

# extracting count data column/header info.
count_df = read_delim(COUNT_FILE, delim = "\t", skip = NSKIP, show_col_types = F)

count_mtrx = count_df  %>% 
            select(grep(".bam", names(count_df)))
colnames(count_mtrx) = str_extract(basename(colnames(count_mtrx)), pattern = "^[A-Za-z0-9]+")
# count_mtrx = as.matrix(count_mtrx)
# rownames(count_mtrx) = count_df$Geneid


# ## This is a function that will remove characters that should not end up in the levels for deseq design
# ## Can Add to this as new ones come up.
# remove_unsafe_characters = function(vec){

#   new_vec = str_replace_all(vec, "%", "p_")
#   return(new_vec)
# }


## Creating colData aka col_map
meta_df = read_tsv(META_TSV, show_col_types = F)
meta_df  %>% view()
experiment_samples = meta_df$experiment_desc
experiment_samples

disease = regmatches(experiment_samples, regexpr("(P001[0-9]*)", experiment_samples, perl = TRUE))
disease <- ifelse(disease == "P001", "Control", "RA")  #all P001 labelled samples are control, those with P001[0-9]+ are rheumatoid arthritis
col_map = tibble("accession" = trimws(meta_df$run_accession), 
                 "sample_id" = trimws(meta_df$sample_accession),
                 "disease" = disease
                )
col_map = col_map  %>% 
    mutate("coef1" = disease,
           "replicate" = ave(seq_along(disease), disease, FUN = seq_along),
           "ref_coef1" = "Control"
           )


## Ensure first column of colData equals the matrix colname
if (any(colnames(count_mtrx) != meta_df$run_accession)){
    cat("Matrix Colnames:", colnames(count_mtrx), "\n")
    cat("ColData Colanmes:", meta_df$study_accession, "\n\n")
    stop("Matrix column names do not align with col_map column1 names")
}

## write col_map to file.
  write.table(
    col_map,    
    file.path(OUT_DIR, "colData.tsv"),
    col.names = TRUE,
    row.names = FALSE
  )

