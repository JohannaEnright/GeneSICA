#!/usr/bin/env python3


### IMPORT LIBRARIES ###
import argparse
import random
import os
import time
from statsmodels.stats.multitest import multipletests
import multiprocessing as mp
import pandas as pd
import numpy as np
import warnings

import pareto_dist
from calculate_SIC_scores import read_in_gmt, read_in_goi, calculate_SIC_full



########################
print(time.ctime())
start_time = time.perf_counter()
random.seed(10)

### GLOBAL VARIABLES FOR PARALELLIZATION ##
# parallelization of permutation value computation
# --- worker-side globals (one copy per worker process) ---
_worker_gene_set_dict = None
_worker_gene_set_db = None

#parallelization of gpd approximation
_worker_sic_perm_dict = None
_worker_transformation_exp = None
###########################################

## functions for permutation paralellization
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
def bootstrap_genesica(gene_set_dict, gene_set_db, sic_scores, m, n, n_jobs, max_tail=500, seed_base=10):

    tasks = [(i, m, seed_base + i) for i in range(n)]
    all_sic_perms = {gs : [0, np.zeros(max_tail+1)] for gs in sic_scores}  #1 row, 501 columns (need 501 for t calc)

    # run bootstrap iterations in parallel
    with mp.Pool(processes=n_jobs, initializer=_init_worker,
                 initargs=(gene_set_dict, gene_set_db)) as pool:
        #_unordered means workers can finish in any order
        for count, result in enumerate(pool.imap_unordered(_bootstrap_task, tasks), 1):
            if count % 100 == 0 or count == n:
                print(f"{count}/{n} bootstraps complete", flush=True)
            ## Add all SIC permutation values to a large dictionary
            
            for gs, score in result.items():
                if score > sic_scores[gs]:
                    all_sic_perms[gs][0] += 1   #computing M value (number of perm values greater than test stat
                arr = all_sic_perms[gs][1]  #numpy array holding up to 501 most extreme values
                if score > np.min(arr):
                    arr[np.argmin(arr)] = score  #replace lowest value with new score

    return all_sic_perms


# This function is fed to the parallel pval calc function
def _calc_pval_task(gs):

    test_stat = _worker_sic_scores[gs]
    transformation_exp = _worker_transformation_exp
    perm_dict = _worker_sic_perm_dict

    #aquire computed gene set score information
    gs_info = perm_dict[gs]
    M = gs_info[0]
    perm_list = gs_info[1]

    n_nonzero = len(perm_list[perm_list != 0])
    follows_gpd = False; c = "NA"; scale = "NA"; N_exc= "NA"; in_range="NA" #iniialize variables

    if M < 10:
        fit_info = pareto_dist_test.determine_Y_exc(
            perm_list**transformation_exp,
            #use most extreme 500 (+1) values, or all non-zero values (perm_list should already be cut to 500
            init_chop=min(n_nonzero, max_tail)
            ) 
        if fit_info:
            follows_gpd = True
            c, t, scale, N_exc = fit_info["c"], fit_info["t"], fit_info["scale"], len(fit_info["y_exc"])

            in_range = pareto_dist_test.check_z_range(c, scale, perm_list, t) #will be outputted to file
            pval = pareto_dist_test.approximate_p_val(
                M, test_stat**transformation_exp, c, t, scale,
                N_perm=n, N_exc=N_exc
            )
        else: #this is if the determine_Y_exc function returns that the data does not follow a gpd function
            pval = pareto_dist_test.approximate_p_val(M, N_perm=n)
    else:
        pval = pareto_dist_test.approximate_p_val(M, N_perm=n)

    return gs, pval, follows_gpd, n_nonzero, c, scale, N_exc, M, in_range



def calculate_pvals_parallel(sic_scores, sic_perm_dict, transformation_exp, n_jobs):

    global _worker_sic_scores, _worker_sic_perm_dict, _worker_transformation_exp
    _worker_sic_scores = sic_scores
    _worker_sic_perm_dict = sic_perm_dict
    _worker_transformation_exp = transformation_exp

    pval_dict = {}
    gs_list = list(sic_perm_dict.keys())
    total = len(gs_list)

    with mp.Pool(processes=n_jobs) as pool:   #global arguments makes it so there is no initializer/initargs
        for count, (gs, pval, follows_gpd, n_nonzero, c, scale, N_exc, M, in_range) in enumerate(
                pool.imap_unordered(_calc_pval_task, gs_list), 1):  #chooses 20 genesets at a time

            try:
                a_div_c = scale/c
            except TypeError:
                a_div_c = "NA"
            pval_dict[gs] = (pval, M, follows_gpd, c, scale, a_div_c, in_range, N_exc) 

            if count % 500 == 0 or count == total:
                print(f"{count}/{total} p-values calculated\n", flush=True)

    return pval_dict


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


#pval dict= gs:(pval, M, follows_gpd, c, scale, scale/c, in_range, N_exc)

def adjust_pvals(p_val_dict, sic_scores, mult_test="fdr_bh"):

    sic_df = pd.DataFrame(sic_scores.items(), columns=["gene_set_name", "sic"]) 
    #pval_df = pd.DataFrame(p_val_dict.items(), columns=["gene_set_name", "pval"]) 
    pval_df = pd.DataFrame({
            "gene_set_name" : p_val_dict.keys(),
            "pval" : [i[0] for i in p_val_dict.values()],
            "M" :    [i[1] for i in p_val_dict.values()],
            "follows_gpd" :[i[2] for i in p_val_dict.values()],
            "k" : [i[3] for i in p_val_dict.values()],
            "a" :[i[4] for i in p_val_dict.values()],
            "a_div_k" : [i[5] for i in p_val_dict.values()],
            "in_range" :[i[6] for i in p_val_dict.values()],
            "N_exc" :[i[7] for i in p_val_dict.values()],
            }) 
    score_df = pd.merge(sic_df, pval_df, how="left", on="gene_set_name")
    score_df = score_df.sort_values(by=["pval", "sic"], ascending=[True,False]) #order by pvalue and sic

    adj_pvals = fdr_bh_calc(score_df["pval"].tolist(), mult_test) #calculate BH adjusted pvals, according to previous ordering
    score_df["p.adj"] = adj_pvals
    #reorder the columns so padj is beside pval
    score_df = score_df[[
        "gene_set_name", "sic", "pval", "p.adj", "M",
        "follows_gpd", "k", "a", "a_div_k", "in_range", "N_exc"
        ]]

    return score_df  #return the whole dataframe to retain the correct order


#new version of write to function which uses the dataframe and not dicitonaries
def write_pval_to_file(score_df, out_dir, gmt_file, goi_file, m, n, transformation_exp, minsize, maxsize):
    
    #basename = f"genesica_pval_sic_order_pseudocount_{m}m_{n}n.tsv"
    basename = f"genesica_pval_gpd_{transformation_exp}exp_{m}m_{n}n_{minsize}min_{maxsize}max.tsv"
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
        "--n_jobs",
        type=int,
        default=1, 
        help="Number of parallel processing jobs"
        )
    parser.add_argument(
        "--max_tail",
        type=int,
        default=500, 
        help="Number of inital permutation data points used to generate the pareto distribution."
        )
    parser.add_argument(
        "--transformation_exp",
        type=float,
        default=1, 
        help="""
        Exponential transformation you want to apply to the test statistic and exceedence values.
        Value must be greater than 0 (only increasing functions allowed).
        Transformations are used to get a heavy-tailed GPD, (shape parameter k < 0)
        Default: 1
        """
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


       
    args = parser.parse_args()

    gmt_file = args.gmt
    goi_file = args.goi_file
    mult_test = args.multiple_test
    out_dir = args.out_dir
    m = args.m
    global n  #making n a global variable for later on.
    n = args.n
    n_jobs = args.n_jobs
    transformation_exp = args.transformation_exp
    global max_tail
    max_tail = args.max_tail
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
    sic_perm_dict = bootstrap_genesica(gene_set_dict, gene_set_db, sic_scores,
                                        m, n, n_jobs=n_jobs, max_tail=max_tail) #gs:[M, np.array(500,1)]

    #parallelization of calculating p-values. gs:(pval, M, follows_gpd, c, scale, scale/c, in_range, N_exc)
    pval_dict = calculate_pvals_parallel(sic_scores, sic_perm_dict, transformation_exp, n_jobs)

    ##########################
    #Calculate p-values after adjusting for multiple hypothesis testing.
    score_df = adjust_pvals(pval_dict, sic_scores, mult_test=mult_test)

    #write p-values for each gene set to a file
    write_pval_to_file(score_df, out_dir, gmt_file, goi_file, m, n, transformation_exp, minsize, maxsize)

    #measure time elapsed
    end_time = time.perf_counter()
    elapsed_time = (end_time - start_time)/3600
    print(f"Elapsed time: {elapsed_time:.4f} hours")


if __name__ == "__main__":
    main()

