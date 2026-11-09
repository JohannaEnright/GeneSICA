#!/usr/bin/env Rscript
#Creating the contrast.txt file

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
parser$add_argument("--count_file",              ## I don't actually do anything with this variable yet.
                    default = "", 
                    type = "character",
                    help = "file path and name of raw expression count matrix. exp. 'expression_datasets/GSE265795/counts/read_counts.txt' "
)
parser$add_argument("--out_dir",              ## I don't actually do anything with this variable yet.
                    default = getwd(), 
                    type = "character",
                    help = "Output directory for normalized count matrix."
)
parser$add_argument("--nskip",              ## I don't actually do anything with this variable yet.
                    default = 1, 
                    type = "integer",
                    help = "Number of lines needed to skip before you get to the actual count matrix. Default is 1"
)
parser$add_argument("--filename_add",              ## I don't actually do anything with this variable yet.
                    default = "", 
                    type = "character",
                    help = "This is if you want to add any text to the file name. Like 'GSE300130' for example."
)
parser$add_argument("--col_map",              ## I don't actually do anything with this variable yet.
                    default = "", 
                    type = "character",
                    help = "This is the colData.tsv file. It contains meta data about the raw count file. (generated from a create_col_data script)"
)
parser$add_argument("--contrast_vars",
                       nargs = 2,
                       type = "character",
                       help = "The 2 variables (relevant to the specific samples/dataset) which you want to be contrasted.  Can be 2 subcategories in 'coef1' for example. ('KIF20AOE' vs 'CTRL')."
)
parser$add_argument("--control_case",
                       nargs = 2,
                       type = "character",
                       help = "Indicate which of the contrast_vars is a control, and which is the case by indicating 'c' for control and 'd' for case. (Order must correspond with contrast_var order)"
)

##### ARGPARSER ####
args = parser$parse_args()
nskip = args$nskip
out_dir = args$out_dir
count_file = args$count_file
filename_add = args$filename_add
col_map = args$col_map
meta_file = args$meta_file
contrast_vars = args$contrast_vars
control_case = args$control_case

# count_file="../expression_datasets/GSE265795/counts/read_counts.txt"
# out_dir="silver/jo_test_data"
# filename_add="GSE265795"
# col_map="../expression_datasets/GSE265795/meta/colData.tsv"
# contrast_vars_1="KIF20AOE"
# contrast_vars_2="CTRL"
# contrast_vars=c("KIF20AOE", "CTRL")
# nskip=1
# control_case=c("d", "c")

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
read_count_file = function(count_file, nskip){

  count_df = read_delim(count_file, delim = "\t", skip = nskip, show_col_types = F)
  
  count_mtrx = count_df  %>% select(grep(".bam", names(count_df)))
  colnames(count_mtrx) = str_extract(basename(colnames(count_mtrx)), pattern = "^[A-Za-z0-9]+") #everything up do \.
  count_mtrx = as.matrix(count_mtrx)
  rownames(count_mtrx) = count_df$Geneid

  return(count_mtrx)
}


## reading in colData.tsv file (meta data for count matrix columns and how design should be specified)
read_in_col_map = function(col_map){

  cat("Loading ColData\n")
  col_map_df = read_delim(col_map, show_col_types = F)

  return(col_map_df)
}


#A subfunction of find_countrast_columns
contrast_rows = function(col_map_df, var){

  coef_cols = col_map_df[grepl("^coef", names(col_map_df))]
  print(coef_cols)
  row_idx = which(apply(coef_cols, 1, function(row) var %in% row))
  filtered_rows = col_map_df[row_idx,]
  return(filtered_rows)

}

## basically, I only want to keep columns in the count matrix that correspond to the test + control samples
##returns a dataframe of the count matrix column I need to keep ('accession') and whether each column is a case or control sample
find_contrast_columns = function(count_mtrx, col_map_df, contrast_vars, control_case){

  c1=contrast_vars[[1]]
  c2=contrast_vars[[2]]
  ## Confirm that both contrast variables are in the meta data
  coef_cols <- grep("^coef", colnames(col_map_df), value = TRUE)

  if (!c1 %in% unlist(col_map_df[, coef_cols])){
    stop(glue("First Contrast Variable {c1} is not in colData.tsv"))
  }
  if (!c2 %in% unlist(col_map_df[, coef_cols])){
    stop(glue("Second Contrast Variable {c2} is not in colData.tsv"))
  }
  #   if (!c1 %in% unlist(col_map_df)){
  #   stop(glue("First Contrast Variable {c1} is not in colData.tsv"))
  # }
  # if (!c2 %in% unlist(col_map_df)){
  #   stop(glue("Second Contrast Variable {c2} is not in colData.tsv"))
  # }

  ## First get the parts of the df that correspond to the contrast variable
  filtered_1 = contrast_rows(col_map_df, c1)
  filtered_2 = contrast_rows(col_map_df, c2)

  #classify which is control and which is case
  filtered_1$contrast = control_case[[1]]
  filtered_2$contrast = control_case[[2]]

  #create a new mapping df.
  map1 = filtered_1  %>% 
           select(accession, contrast)  %>% 
           mutate("var"=c1)
  map2 = filtered_2  %>% 
            select(accession, contrast) %>% 
           mutate("var"=c2)
  contrast_map = bind_rows(map1, map2)

  #determine which column numbers the accession names correspond to in the count matrix
  count_col_idx = match(contrast_map$accession, colnames(count_mtrx))
  contrast_map$col_num = count_col_idx

  #now I will order by col_num, so they are in the same order as the count_matrix.
  contrast_map = contrast_map  %>% arrange(col_num)

  return(contrast_map)

}

write_contrast_map = function(contrast_map, out_dir, filename_add){
  
  if (filename_add != "") {
    file_basename = paste0(filename_add, "_contrast_accession_map.txt")
  } else {
    file_basename = "contrast_accession_map.txt"  
  }
    filename = file.path(out_dir, file_basename)
    write.table(contrast_map, file=filename, sep="\t", quote=F, col.names=TRUE, row.names=FALSE)
    print(glue("Contrast map written to:              {filename}"))

}

write_contrast_file = function(contrast_map, out_dir, filename_add){

  if (filename_add != "") {
    file_basename = paste0(filename_add, "_contrast.txt")
  } else {
    file_basename = "contrast.txt"  
  }
    filename = file.path(out_dir, file_basename)
    cat(paste(contrast_map$contrast, collapse = "\t"), file = filename) ## 'c's and 'd's are in same order as columns in count matrix
    print(glue("Contrast file for silver written to:  {filename}"))
}

###### Main Function --> to call all the others ###
main_func = function(){

  create_output_dir(out_dir)
  count_mtrx = read_count_file(count_file, nskip)
  col_map_df = read_in_col_map(col_map)
  contrast_map = find_contrast_columns(count_mtrx, col_map_df, contrast_vars, control_case)
  write_contrast_map(contrast_map, out_dir, filename_add)
  write_contrast_file(contrast_map, out_dir, filename_add)

}


if (!interactive()) {
  main_func()
  print("FINISHED")
}

warnings()

