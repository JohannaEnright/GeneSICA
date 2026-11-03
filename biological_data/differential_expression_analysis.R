#!/usr/bin/env Rscript

####### LOAD LIBRARIES #########
#BiocManager::install("apeglm")
print("Loading Required Packages...")
suppressPackageStartupMessages({
library (DESeq2) #load the package
library(readr)
library(tidyverse)
library(stringr)
library(argparse)
library(glue)
})
################################
options(warn = 1) #print warning messages as they occur



### Handle creation of output directory 
create_output_dir = function(out_dir){

  if (out_dir == ""){
    out_dir = "out"
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  return(out_dir)
}


### This function reads in count matrix created from .bam files
read_count_file_de = function(count_file, nskip, out_dir){

  count_df = read_delim(count_file, delim = "\t", skip = nskip, show_col_types = F)
  
  count_mtrx = count_df  %>% select(grep(".bam", names(count_df)))
  colnames(count_mtrx) = str_extract(basename(colnames(count_mtrx)), pattern = "^[A-Za-z0-9]+") #everything up do \.
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

  #print(count_mtrx)
  return(count_mtrx)
}


## reference level columns must be named "ref_coef[0-9]" where coef[0-9] corresponds to the correct coefficient column
relevel_factored_coefs = function(col_map_df, coef, ref_level){

  if (length(ref_level) != 1){
    stop(glue("Can only specify 1 reference level, not {length(ref_level)} for {coef}."))
  }

  col_map_df = col_map_df  %>% mutate(
    !!coef := relevel(.data[[coef]], ref_level)
  )
  return(col_map_df)

}

## reading in colData.tsv file (meta data for count matrix columns and how design should be specified)
## if you want to set a reference level comparator, need to specifiy a 'ref_coef[0-9]+' column in colData.tsv
read_in_col_map = function(col_map){

  cat("Loading ColData\n")
  col_map_df = read_delim(col_map, show_col_types = F)

  #adjust the coef columns to be factored with a designated reference level
  col_map_df = col_map_df  %>% mutate(across(starts_with("coef"), ~factor(.x, levels = sort(unique(.x))))) #convert coefficients to factors
  coef_cols = col_map_df  %>% select(starts_with("coef"))
  ref_cols = col_map_df  %>% select(starts_with("ref_coef"))

  if (ncol(ref_cols) == 0){
    warning("No reference levels set. Factoring alphabetically...")
  } else {
    for (i in 1:ncol(coef_cols)){
      coef = names(coef_cols)[i]
      ref_col = paste0("ref_", coef)
      if (!ref_col %in% names(ref_cols)){
        warning(glue("{names(coef)} does not have an assigned reference level.    Factoring alphabetically..."))
        next
      } else {
        ref_level = unique(ref_cols[[ref_col]])
        col_map_df = relevel_factored_coefs(col_map_df, coef, ref_level) #Relevel according to the refernece level
      }
    }
  }
  cat("\nFactor Levels:\n")
  print(col_map_df %>% select(starts_with("coef")) %>% lapply(levels))
  cat("\n")

  return(col_map_df)
}

# use this to determine the interaction order for the design matrix based on the number of replicates in the sample
# it assumes replicates are labelled in the column 'replicate' and coeficients are in 'coef[1-9]+' column
# interaction order will only go as high as ^3 and is determined by the number of replicates for each biological state.
get_interaction_order = function(col_map_df){
  
  cat("Getting interaction order for design...\n")
  coefs_df = dplyr::select(col_map_df, grep("^coef", ignore.case = TRUE, names(col_map_df))) 
  n_coefs = ncol(coefs_df)

  n_levels <- apply(coefs_df, 
                    MARGIN = 2, 
                    FUN = function(x) length(unique(x)))
  min_reps = col_map_df %>% count(across(starts_with("coef"))) %>% pull(n) %>% min()
  total_samples = nrow(col_map_df)
  
  # Number of parameters for a given interaction order
  n_params = function(order) {
    1 + sum(sapply(1:order, function(k) choose(n_coefs, k)))
  }
  
  cat(sprintf("n coefficients: %d\n", n_coefs))
  #cat(sprintf("Total samples: %d\n\n", total_samples))
  
  # Only check up to order 2 maximum (^2)
  max_order <- min(n_coefs, 2)
  feasible_order <- 1
  
  for (order in 1:max_order) {
    params <- n_params(order)
    min_needed <- params * min_reps
    feasible <- total_samples >= min_needed
    
    if (feasible) feasible_order <- order
  }
  
  cat(sprintf("\nRecommended interaction order: %d\n", feasible_order))
  return(feasible_order)
}

## This is to drop interactions where not all subcategories are present (leads to 0s in matrix)
## Used to determine which variables will have interactions in the 'design'
test_complete_interactions = function(col_map_df, test_pairs){

  count_n_reps = col_map_df %>%
    count(across(all_of(test_pairs))) %>%
    complete(!!!syms(test_pairs), fill = list(n = 0))

  if (0 %in% count_n_reps$n){
    return(FALSE)  #This combination of pairs has some missing combinations (cannot use as interaction)
  } else {
    return(TRUE)   #This combination of pairs has no missing combos, can use as interaction
  }
  
}


## Design created from parsing col_map file for all columns beginning with coef[1-9]+.
## design_type can be 'main_effects', 'pairwise', 'all'  (where all will do all interaction terms if there are > 2 coefficients)
create_design = function(col_map_df){
  
  int_order = get_interaction_order(col_map_df)  #gets the interaction order for the design based on number of coefficients and replicates.
  coefs = names(col_map_df)[grep("^coef[0-9]+", names(col_map_df), ignore.case = T)]  #starts with coef

  ## First order intearction
  if (int_order == 1){
    design = paste("~(", paste(coefs, collapse = " + "), ")")

  ## Second Order Interaction
  } else if (int_order == 2){
    coefs_for_design = c(coefs)   #start with 1st order terms
    all_combos = combn(coefs, 2, simplify = FALSE)  #gives all pairwise combinations of coef_cols

    for (i in 1:length(all_combos)){                 #need to check that there are no missing combinations
      no_zeros = test_complete_interactions(col_map_df, all_combos[[i]])
      print(glue("This combo has some missing combinations:  {all_combos[[i]]}.  This combo will not be assessed"))
      if (no_zeros){
        coefs_for_design = c(coefs_for_design, paste(all_combos[[i]], collapse = "*"))
        }
    }
    design = paste("~(", paste(coefs_for_design, collapse = "+"), ")")
  }
  return(as.formula(design))
}

### This is where the differential expression analysis will be performed, list of all result objects will be returned
## Each result element in the list is one pairwise comparison for each coefficient in the GLM
# shrinkage is applied to not overestimate fold changes for lowly expressed genes
perform_de_analysis = function(col_map, design, count_mtrx, alpha = 0.05){

  cat("ColData Headers:\t", colnames(col_map), "\n")
  cat("Design:\t", as.character(design), "\n\n")

  ### Safety check ##
  if (!identical(colnames(count_mtrx), pull(col_map, "accession"))){
    cat("Colum names of count matrix need to be matched to the rownames of colmap.\n")
    cat("ColNames matrix: ", colnames(count_mtrx), "\n\n")
    cat("First Col of col_map: ", pull(col_map, "accession"), "\n")
    stop()
    }

  dds = DESeqDataSetFromMatrix(countData = count_mtrx,  colData = col_map, design = design)

  # Collapse if there are multiple technical replicates for the same specimen
  if (nrow(col_map) > length(unique(col_map$sample_id))){
    #warning("There are multiple technical replicates per specimen, summing replicates for each specimen.")
    warning("There are multiple technical replicates per specimen, consider collapsing replicates.")
    #dds = collapseReplicates(dds, groupby = col_map$sample_id)  #collapse 
  }

  #differential expression analysis
  dds = DESeq(dds)
  res_names = resultsNames(dds)  #this is a list of the coefficients for the linear model (specified in design)

  all_res = list()
  for (n in res_names[2:length(res_names)]){
    print(glue("Getting Results for {n}"))
    #res = lfcShrink(dds, coef = n, type = "apeglm", returnList = T)  #apply shrinkage for lowly expressed genes
    res = results(dds, name = n, alpha = alpha)  #2-tailed p-values
    all_res[[n]] = res
  }

  return(all_res)
}


## extract relevant data from results for ORA and GSEA
gather_gsea_ora_data = function(res, alpha){

  gsea_ranked = as_tibble(res, rownames = "gene") %>% 
    filter(!is.na(`stat`)) %>%
    arrange(desc(`stat`)) %>%    #order by ranked test statistic (keep all genes)
    dplyr::select("gene", "stat")

  ora_sig <- subset(res, padj < alpha) %>%  #Only keep significant genes for ORA
    as_tibble(rownames = "gene") %>% 
    filter(!is.na(padj)) %>%
    dplyr::select(gene)

  log2_fc = as_tibble(res, rownames = "gene") %>% 
    filter(!is.na(log2FoldChange)) %>%
    dplyr::select("gene", "log2FoldChange")
  
  return(list("gsea_ranked"=gsea_ranked, "ora_sig"=ora_sig, "log2_fc"=log2_fc))
}


write_gsa_files = function(out_dir, filename, results, gsea_ranked, ora_sig, log2_fc){

  ### GSEA FILE (2 columns needed)
  gsea_file = file.path(out_dir, paste0("GSEA_", filename, "_genes_of_interest.rnk"))
  write_delim(
    gsea_ranked,
    file = gsea_file,
    delim = "\t",
    col_names = F
  )

  ### ORA ##  
  ##### Input for ORA ##  --> #Just a list of genes of interest####
  ora_goi_file = file.path(out_dir, paste0("ORA_", filename, "_genes_of_interest.txt"))
  write.table(
    ora_sig,    #A list of all the significantly differentiated genes.
    ora_goi_file,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  ## Also need to create a reference gene file for ORA ##
  ora_reference_file = file.path(out_dir, paste0("ORA_reference_genes.txt"))
  write.table(
    rownames(results),  #all genes
    ora_reference_file,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  ## Log 2 Fold Change file for GAGE ##
  gage_log2fc_file = file.path(out_dir, paste0("GAGE_", filename, "_log2FC.txt"))
  write_delim(
    log2_fc,
    file = gage_log2fc_file,
    delim = "\t",
    col_names = F
  )

  ## Maybe I should also write a file with the parameters used.??
  
  print(glue("Created File Output: {gsea_file}"))
  print(glue("Created File Output: {gage_log2fc_file}"))
  print(glue("Created File Output: {ora_goi_file}"))
  print(glue("Created File Output: {ora_reference_file}"))
}


#write all extra data to a file
write_meta_data = function(results, col_map, filename, out_dir){
  
  meta_out = file.path(out_dir, "extra_data")
  dir.create(meta_out, recursive=TRUE, showWarnings=FALSE)

  # Write raw results and meta data to file
  results_with_rownames <- cbind("Gene" = rownames(results), results)
  full_results = file.path(meta_out, paste0("full_results_", filename, ".txt"))

  write.table(
    results_with_rownames,
    full_results,
    col.names = TRUE,
    row.names = FALSE
  )
  result_params = file.path(meta_out, paste0("result_params_", filename, ".txt"))
  write.table(
    mcols(results, use.names = TRUE) ,
    result_params,
    col.names = NA
  )
  print(glue("Created File Output: {full_results}"))
  print(glue("Created File Output: {result_params}"))
}

## Collect appropriate info and write to files.
parse_results = function(all_res, col_map, out_dir, alpha){

  for (r in names(all_res)){
    output = gather_gsea_ora_data(all_res[[r]], alpha)  #function from sourced file
    gsea_ranked = output$gsea_ranked
    ora_sig = output$ora_sig
    log2_fc = output$log2_fc
    rm(output)
    
    write_gsa_files(out_dir, r, all_res[[r]], gsea_ranked, ora_sig, log2_fc)  #write files for gsa in webgestalt
    write_meta_data(all_res[[r]], col_map, r, out_dir)
  }
}


## The head cheif function --> calling all functions
de_main_func = function(){
  
  ############# DIRECTORIES AND FILES --> GLOBAL VARIABLES #######
parser = ArgumentParser(description = "Differential Expression Analysis for RNAseq Data ")

# optional arguments with defaults
parser$add_argument("--count_file",
                    default = "",
                    type = "character",
                    help = "File countaining raw RNA expression counts."  #counts should not be normalized.
)
parser$add_argument("--out_dir",
                    default = "",   #default output is where the count_file is located
                    type = "character",
                    help = "Output Directory for Files"
)
parser$add_argument("--nskip",
                    default = 1,     #skip the first line which contains a note, then you can read in matrix as a delimiter file.
                    type = "integer",
                    help = "Number of rows in the matrix count file to skip before getting to the column heading row. Only necessary if count file is provided"
)
parser$add_argument("--alpha",
                    default = 0.05,
                    type = "numeric",
                    help = "FDR Threshold on adjusted p-value"
)
parser$add_argument("--col_map",
                    default = "",   #
                    type = "character",
                    help = "file containing colData for DEseq. Must contain headers 'sampleid', 'coef[1-9]+', and 'replicate'. Can contain other informatory columns as well."
)


args = parser$parse_args()

nskip = args$nskip
out_dir = args$out_dir
col_map = args$col_map
count_file = args$count_file
alpha = args$alpha

########################################


  out_dir = create_output_dir(out_dir)
  
  ### If input type is matrix of counted read mappings
  count_mtrx = read_count_file_de(count_file, nskip, out_dir)
    
  if (any(count_mtrx %% 1 != 0)){
    warning("Matrix countains fractional amounts: Rounding...")
    count_mtrx = round(count_mtrx, 0)  #round fractional counts (from fractional multimapping reads.)
  }

  #Remove any rows where GeneIds are NA
  before_NA_removal = nrow(count_mtrx)
  count_mtrx = count_mtrx[!is.na(rownames(count_mtrx)), ]
  cat("\n", before_NA_removal-nrow(count_mtrx), " Rows containing GeneIDs as NA were dropped.\n")
  cat(nrow(count_mtrx), "Rows Remaining in Count Matrix\n\n")
  
  
  #Read in colData info
  col_map_df = read_in_col_map(col_map)
  
  ## Determine design of the linear model
  design = create_design(col_map_df)
  
  all_res = perform_de_analysis(count_mtrx = count_mtrx, col_map = col_map_df, design = design, alpha = alpha)

  # Extract valuable information from results
  parse_results(all_res, col_map_df, out_dir, alpha)
}

if (!exists("being_sourced") || !being_sourced) {
  de_main_func()
  warnings()
  print("FINISHED")
}

