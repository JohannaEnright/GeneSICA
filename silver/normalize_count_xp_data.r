#!/usr/bin/env Rscript
# take a gene expression count data matrix and normalizes it. --> This normalization is required for input into Silver

####### LOAD LIBRARIES #########
#BiocManager::install("apeglm")
print("Loading Required Packages...")
suppressPackageStartupMessages({
library(tidyverse)
library(argparse)
library(glue)
library(DESeq2)
})

################################
options(warn = 1) #print warning messages as they occur

############# DIRECTORIES AND FILES --> GLOBAL VARIABLES #######
parser = ArgumentParser(description = "Converting Raw Count Matrix to Normalized Count Matrix")

# optional arguments with defaults
parser$add_argument("--count_file",        
                    default = "", 
                    type = "character",
                    help = "file path and name of raw expression count matrix. exp. 'expression_datasets/GSE265795/counts/read_counts.txt' "
)
parser$add_argument("--out_dir",          
                    default = getwd(), 
                    type = "character",
                    help = "Optional, otherwise will be current directory. Output directory for normalized count matrix."
)
parser$add_argument("--nskip",            
                    default = 1, 
                    type = "integer",
                    help = "Number of lines needed to skip before you get to the actual count matrix. Default is 1"
)
parser$add_argument("--filename_add",           
                    default = "", 
                    type = "character",
                    help = "(Optional). This is if you want to add any text to the file name. Like 'GSE300130' for example."
)
# parser$add_argument("--col_map",              
#                     default = "", 
#                     type = "character",
#                     help = "This is the colData.tsv file. It contains meta data about the raw count file. (generated from a create_col_data script)"
# )
## created from create_contrast_file.r  (It is necessary to filter the normalized count matrix down to just the case and control sample columns ('d' and 'c') required for silver.
parser$add_argument("--contrast_map_file",       
                    default = "", 
                    type = "character",
                    help = "This is the contrast_accession_map.txt. It contains info about which column numbers and names correspond to the count matrix column. It is required for filtering down the count matrix to only relevent case and control variables"
)

args = parser$parse_args()
nskip = args$nskip
out_dir = args$out_dir
count_file = args$count_file
filename_add = args$filename_add
#col_map = args$col_map
contrast_map_file = args$contrast_map_file

# nskip=1
# count_file="../expression_datasets/GSE265795/counts/read_counts.txt"
# out_dir="silver/jo_test_data"
# filename_add="GSE265795"
# col_map="../expression_datasets/GSE265795/meta/colData.tsv"
# contrast_map_file="silver/jo_test_data/GSE265795_contrast_accession_map.txt"



### Create output dir if it doesn't exist
### Handle creation of output directory 
create_output_dir = function(out_dir){

  if (out_dir != ""){
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  } else {
    dir.create("out", recursive = TRUE, showWarnings = FALSE)
  }
  return(out_dir)
}


### Read in count file matrix (created from featureCounts, from .bam files)
read_count_file = function(count_file, contrast_map, nskip, out_dir){

  count_df = read_delim(count_file, delim = "\t", skip = nskip, show_col_types = F)
  
  count_mtrx = count_df  %>% select(grep(".bam", names(count_df)))
  colnames(count_mtrx) = str_extract(basename(colnames(count_mtrx)), pattern = "^[A-Za-z0-9]+") #everything up do \.

  #Select only columns in the contrast map file
  count_mtrx = count_mtrx  %>% 
                  select(contrast_map$accession)
  cat("Selecting columns with accessions:\n")
  cat(paste(contrast_map$accession, collapse="\t"), "\n")
  cat("These contain the case and control variables:\n")
  cat(paste(unique(contrast_map$var), collapse="\t"), "\n")

  #Make it an actual matrix
  count_mtrx = as.matrix(count_mtrx)
  rownames(count_mtrx) = count_df$Geneid


  ## Need to remove rows where Genes come from alternate assemblies/haplotypes and have an _1 or _2 for example
  alternate_haplotypes_ind = grepl("_[0-9]+$", rownames(count_mtrx))
  haps_removed = rownames(count_mtrx[alternate_haplotypes_ind, ])  %>% sort()
  count_mtrx = count_mtrx[!alternate_haplotypes_ind,]

  #write removed gene haplotypes to file.
  alt_haps_filename = file.path(out_dir, "removed_alt_haps_assmbl.txt")
  writeLines(haps_removed, con = alt_haps_filename)
  cat(glue("Removing {length(haps_removed)} alternate haplotypes/assemblies."), "\n")
  cat(glue("Wrote removed haplotypes and assemblies file to {alt_haps_filename}\n"))

  ## Round any fractional counts to the neares whole number
  if (any(count_mtrx %% 1 != 0)){
        warning("Matrix countains fractional amounts: Rounding...")  #this is only if values are from fractional counts of multimapped reads
        count_mtrx = round(count_mtrx, 0)  #round fractional counts (from fractional multimapping reads.)
    }

    #Remove any rows where GeneIds are NA
    before_NA_removal = nrow(count_mtrx)
    count_mtrx = count_mtrx[!is.na(rownames(count_mtrx)), ]
    cat("\n", before_NA_removal-nrow(count_mtrx), " Rows containing GeneIDs as NA were dropped.\n")
    cat(nrow(count_mtrx), "Rows Remaining in Count Matrix\n\n")


  return(count_mtrx)
}



## read in contrast map
read_in_contrast_map = function(contrast_map_file){

  contrast_map_df = read_delim(contrast_map_file, delim="\t", show_col_types = FALSE)
  control = contrast_map_df[contrast_map_df$contrast == "c",]  %>% 
      pull(var)  %>% 
      unique()
  contrast_map_df$var = as.factor(contrast_map_df$var)
  contrast_map_df$var = relevel(contrast_map_df$var, control)

  return(contrast_map_df)

}

#write matrix to a file
write_normal_matrix = function(normalized_counts, out_dir, filename_add){

  if (filename_add != "") {
    file_basename = paste0(filename_add, "_profile.txt")
  } else {
    file_basename = "profile.txt"  
    }
  filename = file.path(out_dir, file_basename)
  Geneid = tibble("ID" = rownames(normalized_counts))
  normalized_counts = as_tibble(normalized_counts)
  normalized_counts = bind_cols(Geneid, normalized_counts)
  write.table(normalized_counts, file=filename, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  print(glue("Normalized Count File Written to:  {filename}"))

}


### Normalize the Count Matrix
normalize_count_matrix = function(count_mtrx, contrast_map_df){

    dds <- DESeqDataSetFromMatrix(countData = count_mtrx, colData = contrast_map_df, design = ~ var)
    dds <- estimateSizeFactors(dds) # Ensure size factors are estimated
    normalized_counts <- log2(counts(dds, normalized = TRUE)+1)  ## Log2 normalization
    return(normalized_counts)
}



###### Main Function --> to call all the others ###
main_func = function(){
  
    out_dir = create_output_dir(out_dir)
    contrast_map_df = read_in_contrast_map(contrast_map_file)
    count_mtrx = read_count_file(count_file, contrast_map_df, nskip, out_dir)
    normalized_counts = normalize_count_matrix(count_mtrx, contrast_map_df)
    write_normal_matrix(normalized_counts, out_dir, filename_add)
}

if (!interactive()) {
  main_func()
  warnings()
  print("FINISHED")
}


