#!/usr/bin/env Rscript
### Differential Expression Analysis for Silver specifically (because it has normalized matrix values, cannot use DESEQ2)


print("Loading Required Packages...")
suppressPackageStartupMessages({
library(tidyverse)
library(stringr)
library(argparse)
library(glue)
library(edgeR)
library(limma)
})
################################
options(warn = 1) #print warning messages as they occur


################# ARGUMENTS ###############

 ############# DIRECTORIES AND FILES --> GLOBAL VARIABLES #######
parser = ArgumentParser(description = "Differential Expression Analysis for RNAseq Data ")

# optional arguments with defaults
parser$add_argument("--count_file",
                    default = "",
                    type = "character",
                    help = "File countaining SILVER simulated expression counts."  
)
parser$add_argument("--out_dir",
                    default = "de_out",
                    type = "character",
                    help = "Output Directory for Files"
)
parser$add_argument("--nskip",
                    default = 0,     #skip the first line which contains a note, then you can read in matrix as a delimiter file.
                    type = "integer",
                    help = "Number of rows in the matrix count file to skip before getting to the column heading row. Only necessary if count file is provided"
)
parser$add_argument("--alpha",
                    default = 0.05,
                    type = "numeric",
                    help = "Critical Value (Default=0.05)"
)
parser$add_argument("--contrast_file",
                    default = "",   #
                    type = "character",
                    help = "'simulated.contrast.txt' outputed by SILVER"
)
args = parser$parse_args()

nskip = args$nskip
out_dir = args$out_dir
contrast_file = args$contrast_file
count_file = args$count_file
alpha = args$alpha

##########################################






### Handle creation of output directory 
create_output_dir = function(out_dir){

  if (out_dir == ""){
    out_dir = "out"
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  return(out_dir)
}


### This function reads in count matrix created from .bam files
read_count_file_desilver = function(count_file, nskip){

  count_df = read_delim(count_file, delim = "\t", skip = nskip, show_col_types = F)
  names(count_df)[1] = "Geneid"
  Geneid = count_df$Geneid
  count_mtrx = as.matrix(count_df  %>% select(-Geneid)) #keep only numerical values in the matrix
  rownames(count_mtrx) = Geneid

  return(count_mtrx)
}

read_in_contrast_file = function(contrast_file){

    contrast_vec = readLines(contrast_file)  %>% strsplit("\t")
    return(contrast_vec[[1]])
}

## Fit linear model to get de results
get_de_results = function(count_mtrx, contrast_vec){

    #create desitgn matrix
    design_mtrx = model.matrix(~0 + contrast_vec)
    colnames(design_mtrx) = substr(colnames(design_mtrx), nchar(colnames(design_mtrx)), nchar(colnames(design_mtrx))) #fix names

    #keep only genes with high enough expression for analysis
    keep = filterByExpr(count_mtrx, design = design_mtrx)  #remove genes with too low expression
    count_mtrx = count_mtrx[keep,]

    #fit initial model
    model = lmFit(count_mtrx, design = design_mtrx)

    #make contrasts
    contrast_matrix <- makeContrasts(d - c, levels = design_mtrx)

    #fit model with contrasts
    model <- contrasts.fit(model, contrast_matrix)
    model = eBayes(model)

    #aquire results
    results = topTable(model, number = Inf)

    return(results)
}



## extract relevant data from results for ORA and GSEA
gather_gsea_ora_data = function(results, alpha=0.05){

  gsea_ranked = as_tibble(results, rownames = "gene") %>% 
    filter(!is.na(t)) %>%
    arrange(desc(t)) %>%    #order by ranked test statistic (keep all genes)
    dplyr::select("gene", "t")

  ora_sig <- subset(results, adj.P.Val <= alpha) %>%  #Only keep significant genes for ORA
    as_tibble(rownames = "gene") %>% 
    filter(!is.na(adj.P.Val)) %>%
    dplyr::select(gene)

  log2_fc = as_tibble(results, rownames = "gene") %>% 
    filter(!is.na(logFC)) %>%
    dplyr::select("gene", "logFC")
  
  return(list("gsea_ranked"=gsea_ranked, "ora_sig"=ora_sig, "log2_fc"=log2_fc))
}


write_gsa_files = function(out_dir, results, gsea_ranked, ora_sig, log2_fc){

  ### GSEA FILE (2 columns needed)
  gsea_file = file.path(out_dir, "GSEA_genes_of_interest.rnk")
  write_delim(
    gsea_ranked,
    file = gsea_file,
    delim = "\t",
    col_names = F
  )

  ### ORA ##  
  ##### Input for ORA ##  --> #Just a list of genes of interest####
  ora_goi_file = file.path(out_dir, "ORA_genes_of_interest.txt")
  write.table(
    ora_sig,    #A list of all the significantly differentiated genes.
    ora_goi_file,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  ## Also need to create a reference gene file for ORA ##
  ora_reference_file = file.path(out_dir, "ORA_reference_genes.txt")
  write.table(
    rownames(results),  #all genes
    ora_reference_file,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  ## Log 2 Fold Change file for GAGE ##
  gage_log2fc_file = file.path(out_dir, "GAGE_log2FC.txt")
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
write_meta_data = function(results, out_dir){
  
  meta_out = file.path(out_dir, "extra_data")
  dir.create(meta_out, recursive=TRUE, showWarnings=FALSE)

  # Write raw results and meta data to file
  results_with_rownames <- cbind("Gene" = rownames(results), results)
  full_results = file.path(meta_out, "full_results.txt")

  write.table(
    results_with_rownames,
    full_results,
    col.names = TRUE,
    row.names = FALSE
  )
  print(glue("Created File Output: {full_results}"))

}


# main function to call other functions
main_func_desilver = function(){

    count_mtrx = read_count_file_desilver(count_file, nskip)
    contrast_vec = read_in_contrast_file(contrast_file)
    results = get_de_results(count_mtrx, contrast_vec)
    result_data = gather_gsea_ora_data(results, alpha)
    out_dir = create_output_dir(out_dir)
    write_gsa_files(out_dir, results,
                result_data$gsea_ranked,
                result_data$ora_sig,
                result_data$log2_fc
                )
    write_meta_data(results, out_dir)
}

main_func_desilver()

