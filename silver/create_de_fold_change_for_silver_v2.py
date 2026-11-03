#!/usr/bin/env python3

## This version 2 script will also create the DE fold change file.
## However, It will take a file with the gene set name \t % upregulated genes in each specific gene set.


### IMPORT LIBRARIES ###
import argparse
import random
import pandas as pd
import json
from pathlib import Path
#######################
random.seed(10)

## This is a list of geneset names which correspond to the names in the GMT file
def read_in_geneset_file(gene_set_file):

    df = pd.read_csv(gene_set_file, sep=r"\s+", engine="python")
    #print(df.head())
    print(df.columns)
   # print(df.shape)
    return df

## Gene set database GMT file
def read_in_gmt(gmt_file):

    gene_set_dict = {}
    with open(gmt_file, "r") as f:
        for line in f.readlines():
            line = line.strip().split("\t")
            genes_in_set = line[2:]
            gene_set_dict[line[0]] = set(genes_in_set)  #create gene set dictionary lookup
    return gene_set_dict


## Need to read in the count profile, so that all DE genes are also present in the profile
def read_in_count_profile(profile_file):

    filtered_genes = []
    with open(profile_file, "r") as f:
        profile = f.readlines()
        for line in profile:
            if "featureCounts" in line:
                continue #filter out
            if "sortedByCoord.out.bam" in line:
                continue #filter out
            gene = line.split("\t")[0]
            filtered_genes.append(gene)

    filtered_genes = set(filtered_genes)
    return filtered_genes


def aquire_genes_in_set(gene_set_name, gene_set_dict, profile_genes):

    if gene_set_name not in gene_set_dict.keys():
            raise NameError(f"Gene set: {gene_set_name}, is not in GMT")
    else:
        gene_set = gene_set_dict[gene_set_name] 
        gene_set = gene_set & profile_genes  #only keep genes found in the profile.
        removed_genes = gene_set - profile_genes  #remove genes if not in the profile.

    return (gene_set, removed_genes)

## Assign a fold change low and high value to each gene.
def assign_FCs(gene_set_df, gene_set_dict, profile_genes):

    FC_dict = {}
    all_removed_genes = []
    #iterate over rows
    for i in range(gene_set_df.shape[0]):
        row = gene_set_df.iloc[i, :]
        genes_in_set, removed_genes = aquire_genes_in_set(row.gene_set_name, gene_set_dict, profile_genes)
        all_removed_genes.extend(removed_genes)

        ## Create sampling list with correct proportion of upregulated and downregulated genes
        n_gois = len(genes_in_set)
        n_up = int(round(n_gois*(row.pct_upregulated/100),0))
        n_down = n_gois - n_up
        sampling_list = ["+"]*n_up + ["-"]*n_down
        random.shuffle(sampling_list)

        #iterate over each gene in set
        for gene in genes_in_set:
            up_or_down = sampling_list.pop()  ### remove and return the last item (sampling without replacement)

            if up_or_down == "+":
                FC_dict[gene] = (int(row.low_FC), int(row.high_FC))
            if up_or_down == "-":
                FC_dict[gene] = (int(-1*row.high_FC), int(-1*row.low_FC))
            #gene sets later in the list with overlapping genes, will replace previous genes. 
            #Thus specified percentages of up and downregulated genes may not be exact in the output

    return (FC_dict, set(all_removed_genes))


def what_gene_sets_do_removed_genes_belong_to(removed_genes, gene_set_dict, gene_set_names):

    removed_gene_dict = {}
    for gene in removed_genes:
        gene_sets_in_removed = []
        for gs_name in gene_set_names:
            if gene in gene_set_dict[gs_name]:
                gene_sets_in_removed.append(gs_name)
        removed_gene_dict[gene] = (len(gene_sets_in_removed), gene_sets_in_removed)
    
    return removed_gene_dict

            
## Write the data to a file.
def write_to_FC_file(FC_dict, out_dir):

    #aquire info for file
    to_file = "ID\tFCLower\tFCUpper\n"
    for k, v in FC_dict.items():
        to_file += k + "\t" + str(v[0]) + "\t" + str(v[1]) + "\n"

    #make out_dir if not existing
    Path(out_dir).mkdir(exist_ok=True, parents=True)

    #write to file
    filename = out_dir + "/DE_gene_fold_change.txt"
    with open(filename, "w") as f:
        f.write(to_file) 
    print(f"Wrote data to file:         {filename}")


## Write data on the actual final up and downregulated proportions for each gene set.
def write_meta_data(gene_set_df, FC_dict, gene_set_dict, profile_genes, removed_genes_dict, out_dir):

    ## Aquiring statistics
    final_gs_stats = {}
    unique_FC_combos = []
    for i in range(gene_set_df.shape[0]):
        row = gene_set_df.iloc[i, :]
        genes_in_set, _ = aquire_genes_in_set(row.gene_set_name, gene_set_dict, profile_genes) #throwaway removed genes
        n_genes = len(genes_in_set)
        num_upregulated = 0
        FC_count = {}
        for gene in genes_in_set:
            gene_data = FC_dict[gene]
            ## Count upregulated genes in gene set
            if gene_data[0] > 0 and gene_data[1] > 0 :
                num_upregulated += 1
            ## Count fold change values in gene set
            if gene_data not in FC_count:
                FC_count[gene_data] = 1
                unique_FC_combos.append(gene_data)
            else:
                FC_count[gene_data] += 1

        num_upregulated = round(num_upregulated/n_genes, 3)
        FC_count = {k: round(v/n_genes, 3) for k, v in FC_count.items()}
        final_gs_stats[row.gene_set_name] = (num_upregulated, FC_count)

    unique_FC_combos = (sorted(list(set(unique_FC_combos))))

    to_write = "gene_set_name\tpct_upregulated\t"
    for combo in unique_FC_combos:
        to_write += str(combo[0]) + "_" + str(combo[1]) + "\t"
    to_write += "\n"
    for gs_name, info in final_gs_stats.items():
        to_write += gs_name + "\t" + str(info[0]) + "\t" 
        for combo in unique_FC_combos:
            if combo in info[1]:
                to_write += str(info[1][combo]) + "\t"
            else:
                to_write += "0\t"
        to_write += "\n"
    
    to_write += "\n\nRemoved Genes (not found in profile):\n" + json.dumps(removed_genes_dict, indent=4, sort_keys=True)


    ## Write to a file 
    filename = out_dir + "/meta_DE_fc_v2.txt"
    with open(filename, "w") as f:
        f.write(to_write) 
    print(f"Wrote data to file:         {filename}")



def main():
    
   #### Use Argparse to Read in Variables ########
    parser = argparse.ArgumentParser(description = 
    """
    This script takes a file of gene set names (column1), percentage of upregulated genes (column 2), low fold change (column 3) and high fold change (column 4). 
    It will return the de_fold_change file with gene differential expression required for input to silver.
    """
    )
    parser.add_argument(
        "--gene_set_file",
        type=str,
        default="silver/jo_test_data/gene_sets_cell_cycle.txt",
        help="file containing gene sets which are wanted to be differentially enriched. And percentage of upregulated genes in that gene set. Followed by the lower bound fold change and upper bound fold change for genes in that gene set. Tab separated file."
        )
    parser.add_argument(
        "--gmt_file",
        type=str,
        default="../gmt/c5.go.v2026.1.Hs.symbols.gmt",
        help="GMT file containing the specified gene sets in 'gene_set_file' "
        )
    parser.add_argument(
        "--profile_file",
        type=str,
        default="", #../expression_datasets/GSE300130/counts/reverse_strand/counts/read_counts.txt 
        help="Expression profile used to ensure genes in DE are present in profile."
        )
    parser.add_argument(
        "--out_dir",
        type=str,
        default="silver/jo_test_data", 
        help="output directory for 'DE_gene_fold_change.txt' file"
        )

    args = parser.parse_args()
    gene_set_file = args.gene_set_file
    gmt_file = args.gmt_file
    out_dir = args.out_dir
    profile_file = args.profile_file
    ############################# End of Argparse ################

    print(f"Gene Set File:              {gene_set_file}")
    print(f"GMT File:                   {gmt_file}")
    gene_set_df = read_in_geneset_file(gene_set_file)  #list of gene set names of interest
    gene_set_dict = read_in_gmt(gmt_file)
    profile_genes = read_in_count_profile(profile_file) #list of genes found in profile (read counts)
    FC_dict, removed_genes = assign_FCs(gene_set_df, gene_set_dict, profile_genes)
    removed_genes_dict = what_gene_sets_do_removed_genes_belong_to(removed_genes, gene_set_dict, gene_set_df.gene_set_name.tolist())

    write_to_FC_file(FC_dict, out_dir)
    write_meta_data(gene_set_df, FC_dict, gene_set_dict, profile_genes, removed_genes_dict, out_dir)


if __name__ == "__main__":
    main()

