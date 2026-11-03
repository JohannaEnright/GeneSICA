#!/usr/bin/env python3

## This script takes a gmt file and a list of genes of interest (1 line per gene, as per ORA protocol)
## It calculates the Shared Information Content for each gene set and returns it as a dictionary {gene_set_name: SIC}
## It does this per the GeneSICA method


######## import libraries ##########
import argparse
import numpy as np
import warnings






####### FUNCTIONS ##### 

#reads in file with 1 gene per line (ORA genes of interest) ('L'). Returns set of genes {gene1, gene2}
def read_in_goi(goi_file):

    with open(goi_file, "r") as f:
        goi_list = {gene.strip() for gene in f.readlines()}

    return goi_list


#read in gmt (gene set database) file. return a dictionary of gene set name : {gene1, gene2, etc}
# It also returns a set of all genes present in the database {gene1, gene2, etc.}
def read_in_gmt(gmt_file, minsize, maxsize):

    gene_set_dict = {}
    gene_set_db = []
    with open(gmt_file, "r") as f:
        for line in f.readlines():
            line = line.strip().split("\t")
            genes_in_set = set(line[2:])
            if (len(genes_in_set) < minsize) or (len(genes_in_set) > maxsize):
                #warnings.warn(f"Skipping Gene set {line[0]}. Has a size of {len(genes_in_set)}\n")
                continue
            gene_set_dict[line[0]] = genes_in_set  #create gene set dictionary lookup
            gene_set_db.extend(genes_in_set)            #create full set of genes in database
    gene_set_db = set(gene_set_db)
    return (gene_set_dict, gene_set_db)




# presence score takes a dictionary of gene sets {gene_set_name: {gene list}, etc.} and a list of genes
# instead of calculating 1 at a time, I do the calculations for all gene sets in the database at each step
def calculate_presence_score(gene_set_dict, gene_list):

    #print("Calculating Presence Scores")
    presence_scores = {}
    for gene in gene_list:
        gene_score = 0
        for gene_set in gene_set_dict.values():
            if gene in gene_set:
                gene_score += 1   #add 1 to the presence score if the gene is in the gene set
        presence_scores[gene] = gene_score

    return presence_scores


# Calculate information content scores
def calculate_info_score(presence_scores, n_gene_sets):

    #print("Calculating Information Scores")
    info_scores = {}
    for gene, score in presence_scores.items():
        if score == 0:
            info_score = 0
        else:
            info_score = np.log(n_gene_sets/score)
        info_scores[gene] = info_score

    return info_scores


#Calculating Shared Information content --> 
def calculate_shared_info_content(info_scores, gene_set_dict, goi_list):

    #print("Calculating Shared Information Content")
    shared_info_scores = {}
    for gs_name, gene_set in gene_set_dict.items():
        
        numerator = 0
        intersect_gs = goi_list & gene_set
        for gene in intersect_gs:
            numerator += info_scores[gene]
        denominator = 0
        union_gs = goi_list | gene_set
        for gene in union_gs:
            denominator += info_scores[gene]

        sic = numerator/denominator
        shared_info_scores[gs_name] = sic

    return shared_info_scores


## Call all the functions for a full SIC calculation
def calculate_SIC_full(goi_list, gene_set_db, gene_set_dict):

    full_gene_set = goi_list | gene_set_db    # union of L and the genes in the database to get a full set of genes to calculate scores for
    #goi_list = goi_list & gene_set_db  #only keep genes that are actually in the gmt file.
    presence_scores = calculate_presence_score(gene_set_dict, full_gene_set)
    info_scores = calculate_info_score(presence_scores, len(gene_set_dict))
    shared_info_scores = calculate_shared_info_content(info_scores, gene_set_dict, goi_list)
    return shared_info_scores


#### GLOBAL VARIABLES ###
def main():
    
    #### Use Argparse to Read in Variables ########
    parser = argparse.ArgumentParser(description = 
    """
    Calculates Shared Information Content per the GeneSICA method. Requires a GMT file and gene of interest text file.
    """
    )
    parser.add_argument(
        "--gmt",
        type=str,
        default="../input_files/original_subset.gmt",  #gmt, gene set database file
        help="GMT file containing gene sets"
        )
    parser.add_argument(
        "--goi_file",
        type=str,
        default="../input_files/ORA_coef2_hydroponics_vs_15p_PEG_6000_refseq_gene_genes_of_interest.txt",  #this is the file containing 'L', the differentially expressed genes of interest
        help="Gene of interest .txt file with each gene on a new line"
        )
       
    args = parser.parse_args()
    gmt_file = args.gmt
    goi_file = args.goi_file

    if gmt_file == "":
        NameError("A GMT File Needs to be Specified")
    if goi_file == "":
        NameError("A gene of interest .txt file Needs to be Specified")
    ##### END OF ARGPARSE, VARIABLES AQUIRED. ###

    #read in files
    goi_list = read_in_goi(goi_file)
    gene_set_dict, gene_set_db = read_in_gmt(gmt_file, 5, 500)     #returns dict and set: ({gene_set_name: {g1, g2, g3, etc.}} , {g1, g2, g3, etc})
    sic_scores = calculate_SIC_full(goi_list, gene_set_db, gene_set_dict)
    print("DONE")


if __name__ == "__main__":
    main()
