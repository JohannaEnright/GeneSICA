#!/usr/bin/env python3



### IMPORT LIBRARIES ###
import argparse
import random
import os
import time
from calculate_SIC_scores import read_in_gmt, read_in_goi, calculate_SIC_full
from statsmodels.stats.multitest import multipletests
import multiprocessing as mp
import pandas as pd



########################
print(time.ctime())
start_time = time.perf_counter()
random.seed(10)

### GLOBAL VARIABLES FOR PARALELLIZATION ##
# --- worker-side globals (one copy per worker process) ---
_worker_gene_set_dict = None
_worker_gene_set_db = None
###########################################

def _init_worker(gsd, gsdb):
    global _worker_gene_set_dict, _worker_gene_set_db
    _worker_gene_set_dict = gsd
    _worker_gene_set_db = gsdb

def _bootstrap_task(task):
    i, m, seed = task
    random.seed(seed)                      # unique seed per bootstrap iteration
    new_L = generate_gene_set(_worker_gene_set_dict, m)
    return calculate_SIC_full(new_L, _worker_gene_set_db, _worker_gene_set_dict)


#takes gene set dictionary and an integer (length of desired generated random gene set list
# here we are generating L (the gene of interest list)
# This corresponds to algorithm 2 in GeneSICA
def generate_gene_set(gene_set_dict, m):

    S = set()
    gs_names = list(gene_set_dict.keys())
    while len(S) < m:
        rand_gs = random.choice(gs_names)
        rand_gene = random.choice(list(gene_set_dict[rand_gs]))
        S.add(rand_gene)
        #print(rand_gs, rand_gene)

    return S
    



## this function calculates SIC scores for randomly generated Ls
## m is the length desired for L'
## gene_set_db is a set {} containing all unique genes from the gene sets in gene_set_dict
## gene_set_dict is a dictionary containing all gene set names and their corresponding gene sets.
## parallel_jobs
## returns dictionary of gene set name : p-value 
def bootstrap_genesica(gene_set_dict, gene_set_db, sic_scores, m, n, n_jobs, seed_base=10):

    tasks = [(i, m, seed_base + i) for i in range(n)]
    counts = {gs: 0 for gs in sic_scores} #initialize counts for all gene sets to be 0

    # run bootstrap iterations in parallel
    with mp.Pool(processes=n_jobs, initializer=_init_worker,
                 initargs=(gene_set_dict, gene_set_db)) as pool:
        #_unordered means workers can finish in any order
        for count, result in enumerate(pool.imap_unordered(_bootstrap_task, tasks), 1):
            if count % 100 == 0 or count == n:
                print(f"{count}/{n} bootstraps complete", flush=True)
            #check if score is less than actual SIC score. (used for p-val calc later)
            for gs, score in result.items():  
                if sic_scores[gs] <= score:
                    counts[gs] += 1

    return {gs: (counts[gs] + 1)/(n + 1) for gs in sic_scores} #p-val calculated  #add pseudo-count to avoid p-values = 0


#pval list must be ordered from most significant to least significant p-value
#the critical value is not the adjusted pvalue, need to compare the pval to crit value to determine if reject or not
#def fdr_bh_calc(pval_list, mult_test, Q=0.05):
#
#    if mult_test == "fdr_bh":
#        m = len(pval_list)
#        adj_pval_list = []
#        for i in range(len(pval_list)):
#            crit_val = Q*i/m
#            adj_pval_list.append(crit_val)
#
#    return adj_pval_list

## updated FDR BH with corrected adjusted p-values
def fdr_bh_calc(pval_list, mult_test):

    if mult_test == "fdr_bh":
        m = len(pval_list)
        adj_pval_list = []
        for i, pval in enumerate(pval_list):
            adj_pval = min(1, pval*m/(i+1))  #cannot be above 1 (rank+1 to account for 0 indexing)
            adj_pval_list.append(adj_pval)
        
        #Enforce monotonicity (iterate back through)
        for i in range(len(adj_pval_list)-2, -1, -1): #skip the first one as none come before it.
            adj_pval_list[i] = min(adj_pval_list[i], adj_pval_list[i+1])

    return adj_pval_list



##calculate p-values for the SIC score for each gene set
##requires original scores and bootstrapped scores
##takes the get the 1000 bootstrapped SIC scores for each gene set.
#def calculate_p_values(sic_scores, all_sics_prime):
#
#    print("\nCalculating p-values.")
#    p_val_dict = {}
#    for gene_set in sic_scores.keys():
#        counter = 0
#        for sic_dict in all_sics_prime:  
#            if sic_scores[gene_set] <= sic_dict[gene_set]:  #if original score is less or equal to bootstraped score
#                counter += 1
#        p_val_dict[gene_set] = counter/len(all_sics_prime) #length should be equal to n (list of 1000 dicts)  
#
#    return p_val_dict 

# This way of adjusting p-values assumes that they are not already sorted (thus ranks tied p-values as the same)
#def adjust_pvals(p_val_dict, mult_test="fdr_bh"):
#
#    #maybe order by p-value and secondarly by SIC score before 
#    all_gs = [gs for gs in p_val_dict.keys()]
#    all_pvals = [p for p in p_val_dict.values()]
#    _, adj_pvals, _, _ = multipletests(all_pvals, method=mult_test) #alpha=0.05 on default (is_sorted = False)
#    return dict(zip(all_gs, adj_pvals))

# This way of adjusting p-values will take into consideration the ordering of SIC values as well.
# Will be ordered by p-value and then also by SIC value
def adjust_pvals(p_val_dict, sic_scores, mult_test="fdr_bh"):

    sic_df = pd.DataFrame(sic_scores.items(), columns=["gene_set_name", "sic"]) 
    pval_df = pd.DataFrame(p_val_dict.items(), columns=["gene_set_name", "pval"]) 
    score_df = pd.merge(sic_df, pval_df, how="left", on="gene_set_name")
    score_df = score_df.sort_values(by=["pval", "sic"], ascending=[True,False]) #order by pvalue and sic

    #_, adj_pvals, _, _ = multipletests(
    #        score_df["pval"].tolist(),
    #        is_sorted=True,
    #        method=mult_test
    #        ) #alpha=0.05 on default (is_sorted = False)
    adj_pvals = fdr_bh_calc(score_df["pval"].tolist(), mult_test) #calculate BH adjusted pvals, according to previous ordering
    score_df["p.adj"] = adj_pvals

    return score_df  #return the whole dataframe to retain the correct order



#def write_pval_to_file(p_val_dict, sic_scores, adj_pvals, out_dir, gmt_file, goi_file, m, n):
#    
#    basename = f"genesica_pval_sic_order_pseudocount_{m}m_{n}n.tsv"
#    filename = os.path.join(out_dir, basename)
#
#    to_write = f"#geneSICA pvalues generated from:\n#{gmt_file}\n#{goi_file}\n#n = {n}\n#m = {m}\ngene_set_name\tsic\tpval\tp.adj\n"
#    for gs in p_val_dict.keys():
#        to_write += gs + "\t" + str(sic_scores[gs]) + "\t" + str(p_val_dict[gs]) + "\t" + str(adj_pvals[gs]) + "\n"
#        
#    with open(filename, "w") as f:
#        f.write(to_write)
#    print(f"p-value data written to: {filename}")

#new version of write to function which uses the dataframe and not dicitonaries
def write_pval_to_file(score_df, out_dir, gmt_file, goi_file, m, n, minsize, maxsize):
    
    #basename = f"genesica_pval_sic_order_pseudocount_{m}m_{n}n.tsv"
    basename = f"genesica_pval_sic_order_pseudocount_{m}m_{n}n_{minsize}min_{maxsize}max.tsv"
    filename = os.path.join(out_dir, basename)

    note = f"#geneSICA pvalues generated from:\n#{gmt_file}\n#{goi_file}\n#n = {n}\n#m = {m}\n"
    with open(filename, "w") as f:
        f.write(note)
        score_df.to_csv(f, sep="\t", index=False)
    print(f"p-value data written to: {filename}")



def main():
    
   #### Use Argparse to Read in Variables ########
    parser = argparse.ArgumentParser(description = 
    """
    GeneSICA gene enrichment method
    """
    )
    parser.add_argument(
        "--gmt",
        type=str,
        default="",     #"original_subset.gmt",  #gmt, gene set database file
        help="GMT file containing gene sets"
        )
    parser.add_argument(
        "--goi_file",
        type=str,
        default="",  #this is the file containing 'L', the differentially expressed genes of interest
        help="Gene of interest .txt file with each gene on a new line"
        )
    parser.add_argument(
        "--multiple_test",
        type=str,
        default="fdr_bh", 
        help="Method for multiple hypothesis testing. Default ['fdr_bh']"
        )
    parser.add_argument(
        "--out_dir",
        type=str,
        default="out/pvalues",  #this is the file containing 'L', the differentially expressed genes of interest
        help="Gene of interest .txt file with each gene on a new line"
        )
    parser.add_argument(
        "-m",
        type=int,
        default=20, 
        help="An integer Number. Cardinality (Length) of Desired Gene Set."
        )
    parser.add_argument(
        "-n",
        type=int,
        default=1000, 
        help="An integer Number. Number of times L' is resampled to calculated p-values."
        )
    parser.add_argument(
        "--minsize",
        type=int,
        default=5, 
        help="Minimum size for a gene set to be considered."
        )
    parser.add_argument(
        "--maxsize",
        type=int,
        default=500, 
        help="Maximum size for a gene set to be considered."
        )
    parser.add_argument(
        "--n_jobs",
        type=int,
        default=1, 
        help="Number of parallel processing jobs"
        )
       
    args = parser.parse_args()

    gmt_file = args.gmt
    goi_file = args.goi_file
    mult_test = args.multiple_test
    out_dir = args.out_dir
    m = args.m
    n = args.n
    n_jobs = args.n_jobs
    maxsize = args.maxsize
    minsize = args.minsize

    if gmt_file == "":
        raise NameError("A GMT File Needs to be Specified")
    if goi_file == "":
        raise NameError("A Gene of Interest .txt File Needs to be Specified")
    if m == "":
        raise NameError("Gene Set Cardinality Must Be Specified")
    os.makedirs(out_dir, exist_ok=True)

    ##### END OF ARGPARSE, VARIABLES AQUIRED. ###

    #read in files
    goi_list = read_in_goi(goi_file)    #{gene1, gene2, ...}  set of genes
    gene_set_dict, gene_set_db = read_in_gmt(gmt_file, minsize, maxsize)   #dict of gene set and set of all genes in dict

    #check that at least some of the genes from the goi list are in the gene_set_db
    n_genes_in_db = len(gene_set_db & goi_list)
    if n_genes_in_db == 0:
        raise ValueError("None of the genes in the gene of interest file are present in the GMT database")
    print(f"{n_genes_in_db}/{len(goi_list)} ({round(100*n_genes_in_db/len(goi_list),1)}%) genes of interest are in the GMT database")
    print(f"The GMT database contains {len(gene_set_db)} unique genes")

    #calculate OG SIC score
    print("Calculations For Original Gene of Interest List...")
    sic_scores = calculate_SIC_full(goi_list, gene_set_db, gene_set_dict)

    #calculate SIC scores for bootstrapped Ls.
    print("\nCalculations For Bootstrapped Gene of Interest Lists...")
    pval_dict = bootstrap_genesica(gene_set_dict, gene_set_db, sic_scores, m, n, n_jobs)

    #compare OG to Bootstrapped and calculated p-values
    #pval_dict = calculate_p_values(sic_scores, all_sics_prime)

    #Calculate p-values after adjusting for multiple hypothesis testing.
    #adj_pvals = adjust_pvals(pval_dict, sic_scores, mult_test=mult_test)
    score_df = adjust_pvals(pval_dict, sic_scores, mult_test=mult_test)

    #write p-values for each gene set to a file
    #write_pval_to_file(pval_dict, sic_scores, adj_pvals, out_dir, gmt_file, goi_file, m, n)
    write_pval_to_file(score_df, out_dir, gmt_file, goi_file, m, n, minsize, maxsize)

    #measure time elapsed
    end_time = time.perf_counter()
    elapsed_time = (end_time - start_time)/3600
    print(f"Elapsed time: {elapsed_time:.4f} hours")


if __name__ == "__main__":
    main()

