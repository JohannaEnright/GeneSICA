#!/usr/bin/env python3


### IMPORT LIBRARIES ###
import argparse
import json
from pathlib import Path
import pandas as pd
import os



# returns gmt file as a dictionary of {gene_set_name:{gene1, gene2}}
def read_in_gmt(gmt_file):

    gene_set_dict = {}
    with open(gmt_file, "r") as f:
        for line in f.readlines():
            line = line.strip().split("\t")
            genes_in_set = line[2:]
            gene_set_dict[line[0]] = set(genes_in_set)  #create gene set dictionary lookup

    return gene_set_dict


def read_in_overlap_score_json(overlap_score_file):

    if not Path(overlap_score_file).is_file():
        raise NameError(f"{overlap_score_file} does not exist")

    print(f"Loading:    {overlap_score_file}")
    with open(overlap_score_file, "r") as f:
        overlap_score_dict = json.load(f)

    return overlap_score_dict


## Takes a gene set dictionary,
def determine_gene_overlap_jaccard(gene_set_dict):

    print("Calculating Overlap Score")
    overlap_scores = {}
    n_gene_sets = len(gene_set_dict)
    #check presence score of each gene of a given gene sets in all other gene sets.
    # Will divide this sum by the total number of gene sets
    count = 0
    min_score = 0
    max_score = 0
    for gs_name, genes1 in gene_set_dict.items():
        summed_jaccard = 0
        for genes2 in gene_set_dict.values():
            numerator = genes1 & genes2
            denominator = genes1 | genes2
            summed_jaccard += len(numerator)/len(denominator)
        #Now add the summed score to the dictionary for each gene set
        normalized_score = round(1000*summed_jaccard/n_gene_sets, 4)
        if count == 0:
            min_score = normalized_score
            max_score = normalized_score
        else:
            if normalized_score < min_score:
                min_score = normalized_score
            if normalized_score > max_score:
                max_score = normalized_score
        overlap_scores[gs_name] = normalized_score
        count += 1
        if count % 1000 == 0:
            print(f"Jaccard calculated for {count}/{n_gene_sets} gene sets")

    print(f"Min overlap score: {min_score}")
    print(f"Max overlap score: {max_score}")
    return overlap_scores

## Filters gene sets based on the size range and overlap score.
## Only return ones within the ranges (as a dictionary)
def filter_gene_sets(gene_set_dict, overlap_scores, size_range, overlap_range):

    to_keep = {}
    for gs_name, genes in gene_set_dict.items():
        score = overlap_scores[gs_name]
        if (score >= overlap_range[0]) and (score <= overlap_range[1]):
            n_genes = len(genes)
            if n_genes >= size_range[0] and n_genes <= size_range[1]:
                to_keep[gs_name] = (n_genes, score)

    return to_keep

#Write overlap scores to file, so file can be used as intermediate step so they are not calculated every time
def write_overlap_genes_to_file(overlap_scores, gmt_file, out_dir):

    os.makedirs(out_dir, exist_ok=True)
    basename_gmt = Path(gmt_file).stem
    filename = f"{out_dir}/overlap_scores_{basename_gmt}.json"

    with open(filename, "w") as f:
        json.dump(overlap_scores, f, indent=4)
    print(f"Overlap score written to:  {filename}")
        

#Wr
def write_filtered_genes_to_file(gs_filtered, gmt_file, out_dir, size_range, overlap_range):

    os.makedirs(out_dir, exist_ok=True)

    # Build a DataFrame directly instead of a string
    df = pd.DataFrame(
        [(k, v[0], v[1]) for k, v in gs_filtered.items()],
        columns=["gene_set_name", "size", "overlap_score"]
    )

    # Sort by overlap_score, lowest to highest
    df = df.sort_values(by=["overlap_score", "size"], ascending=True)
    #to_write = "gene_set_name\tsize\toverlap_score\n"
    #for k, v in gs_filtered.items():
    #    to_write += k + "\t" + str(v[0]) + "\t" + str(v[1]) + "\n"
    
    basename_gmt = Path(gmt_file).stem
    filename = f"{out_dir}/filtered_{basename_gmt}_{size_range[0]}_{size_range[1]}size_{overlap_range[0]}_{overlap_range[1]}overlap.tsv"

    df.to_csv(filename, sep="\t", index=False)
#    with open(filename, "w") as f:
#        f.write(to_write)
    print(f"Filtered gene set file written to:  {filename}")
        

def main():

    parser = argparse.ArgumentParser(description = 
    """
    This script returns gene sets within a certain size range. 
    """
    )
    parser.add_argument(
        "--gmt",
        type=str,
        default="c5.go.v2026.1.Hs.symbols.gmt", 
        help="GMT file containing gene sets"
        )
    parser.add_argument(
        "--size_range",
        type=int,
        nargs = 2,
        help="Range of gene set sizes you want returned (exp. 4 20 would return gene sets with sizes >=4 and <=20)"
        )
    parser.add_argument(
        "--overlap_range",
        type=float,
        nargs = 2,
        help="Range of gene sets of specified overlaps you want returned (exp. 0.5 2 would return gene sets with overlap scores >=0.5 and <=2)"
        )
    parser.add_argument(
        "--overlap_score_file",
        type=str,
        default="",
        help="File containing pre-calculated overlap scores. If not present, overlap scores will be calculated by default"
        )
    parser.add_argument(
        "--out_dir",
        type=str,
        default="gene_sets_filtered",
        help="output_directory of file"
        )




    args = parser.parse_args()
    gmt_file = args.gmt
    size_range = args.size_range
    overlap_range = args.overlap_range
    overlap_score_file = args.overlap_score_file
    out_dir = args.out_dir
    
    print(f"GMT file:      {gmt_file}")
    print(f"Size range:    {size_range}")
    print(f"Overlap range: {overlap_range}")
    
    gene_set_dict = read_in_gmt(gmt_file)
    if overlap_score_file == "":
        overlap_score_dict = determine_gene_overlap_jaccard(gene_set_dict)
        write_overlap_genes_to_file(overlap_score_dict, gmt_file, out_dir)
    else:
        #check file exists
        overlap_score_dict = read_in_overlap_score_json(overlap_score_file)

    gs_filtered = filter_gene_sets(gene_set_dict, overlap_score_dict, size_range, overlap_range)
    write_filtered_genes_to_file(gs_filtered, gmt_file, out_dir, size_range, overlap_range)





if __name__ == "__main__":
    main()


