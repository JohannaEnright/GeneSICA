"""This module shows an example for simulating expression profiles.

"""
from silver.utils import simulate
import time
import argparse

#Time
print(time.ctime())
start_time = time.perf_counter()

########## ARGUMENTS #######
parser = argparse.ArgumentParser(description = 
    """
    This file runs silver to simulate a gene expression profile.
    """
    )
parser.add_argument(
    "--profile",
    type=str,
    default="", 
    help="Normalized gene expression profile 'profile.txt' "
    )
parser.add_argument(
    "--contrast",
    type=str,
    default="", 
    help="Path to contrast file: 'contrast.txt' "
    )
parser.add_argument(
    "--gene_fc",
    type=str,
    default="", 
    help="Path to fold change file : 'DE_gene_fold_change.txt' "
    )
parser.add_argument(
    "--lower_bound_col_name",
    type=str,
    default="FCLower", 
    help="Lower bound of column name in 'DE_gene_fold_change.txt' "
    )
parser.add_argument(
    "--upper_bound_col_name",
    type=str,
    default="FCUpper", 
    help="Upper bound of column name in 'DE_gene_fold_change.txt' "
    )
parser.add_argument(
    "--fc_id_col_name",
    type=str,
    default="ID", 
    help="Column name of gene id in 'DE_gene_fold_change.txt' "
    )
parser.add_argument(
    "--profile_id_col_name",
    type=str,
    default="ID", 
    help="column name of genes in 'profile.txt.' "
    )
parser.add_argument(
    "--n_simulated_ctrls",
    type=int,
    default=20, 
    help="Number of desired simulated control samples"
    )
parser.add_argument(
    "--n_simulated_cases",
    type=int,
    default=20, 
    help="Number of desired simulated case samples"
    )
parser.add_argument(
    "--out_dir",
    type=str,
    default="./", 
    help="Output directory path"
    )   
parser.add_argument(
    "--replacement",
    type=str,
    default="False", 
    help="Sampling with Replacement [True or False]"
    )   


args = parser.parse_args()

# Input informaiton
PROFILE_ADDRESS = args.profile
CONTRAST_ADDRESS = args.contrast
GENE_FC_ADDRESS = args.gene_fc

# Input column Names
LOWER_BOUND_COL_NAME = args.lower_bound_col_name
UPPER_BOUND_COL_NAME = args.upper_bound_col_name
FC_ID_COL_NAME = args.fc_id_col_name
PROFILE_ID_COL_NAME = args.profile_id_col_name

# Simulation parameters
NUM_SIMULATED_CTRLS = args.n_simulated_ctrls
NUM_SIMULATED_CASES = args.n_simulated_cases
REPLACEMENT = args.replacement
True if REPLACEMENT in ["TRUE", "True", "T"] else False #convert string to boolean

# Output info
OUT_DIR = args.out_dir

# Simulate control and case samples
print("Simulating Expression Profile")
sim_ctrls, sim_cases = simulate(profile_address=PROFILE_ADDRESS,
                                contrast_address=CONTRAST_ADDRESS,
                                gene_fc_address=GENE_FC_ADDRESS,
                                fc_id_col_name=FC_ID_COL_NAME,
                                lower_bound_col_name=LOWER_BOUND_COL_NAME,
                                upper_bound_col_name=UPPER_BOUND_COL_NAME,
                                num_simulated_ctrls=NUM_SIMULATED_CTRLS,
                                num_simulated_cases=NUM_SIMULATED_CASES,
                                profile_id_col_name=PROFILE_ID_COL_NAME,
                                num_repository_reps=10,
                                ctrl_symbol='c',
                                case_symbol='d',
                                profile_sep='\t',
                                contrast_sep='\t',
                                fold_change_sep='\t',
                                alpha=0.05,
                                random_state=123456,
                                sampling_with_replacement=REPLACEMENT
                                )
# Combine the simulated controls and cases
expression_dataset = sim_ctrls.concat(sim_cases)
# Get the expression profile
result = expression_dataset.profile
# Write the simulated dataset to a file
simulated_profile_address = OUT_DIR + '/simulated.profile.txt'
result.to_csv(simulated_profile_address, sep='\t')
print(f"Wrote simulated profile to: {simulated_profile_address}")

# Make a contrast for the simulated dataset and write it to a file
simulated_contrast_address = OUT_DIR + '/simulated.contrast.txt'
with open(simulated_contrast_address, 'w') as fout:
    fout.write('{}\n'.format('\t'.join(['c'] * NUM_SIMULATED_CTRLS +
                                       ['d'] * NUM_SIMULATED_CASES)))
print(f"Wrote contrast file to: {simulated_contrast_address}")


#measure time elapsed
end_time = time.perf_counter()
elapsed_time = (end_time - start_time)/360
print(f"Elapsed time: {elapsed_time:.4f} hours")
