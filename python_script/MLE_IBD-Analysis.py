import numpy as np
import pandas as pd
import re
from pathlib import Path
import matplotlib.pyplot as plt
import socket as socket
import os as os
import sys as sys
import itertools as it
import multiprocessing as mp
import argparse 
sys.path.append('/home/eychenne@newcefe.newage.fr/Documents/BIOINFO/IBD-Analysis/ibd_analysis/analysis_popres/')
from mle_multi_run import MLE_analyse
from mle_estim_error import MLE_estim_error
import sys
from contextlib import redirect_stdout

### Define arguments ###
# Initialize parser
parser = argparse.ArgumentParser(description="Run MLE analysis on IBD data")

# Add arguments with a default value or not 
parser.add_argument("--input_file", "-i", type=str, required=True, help="Path to the IBD-Analysis IBD file")
parser.add_argument("--output_folder", "-o", type=str, required=True, help="Path for output ")
parser.add_argument("--loc_file", "-loc", type=str, required=True, help="Path to the IBD-Analysis coordinates file")
parser.add_argument("--start_params", "-p", type=float, nargs=2, default=[1.0, 1.0], help="Start parameters for MLE, given as a string. Example: 1.0, 1.0")
parser.add_argument("--g", "-g", type=float, default=100.0, help="Value of genome size")
parser.add_argument("--model", "-m",type=str, default="constant", help="Model type (default:constant)")
parser.add_argument("--min_len", type=float, default=2.0, help="Minimum IBD length")
parser.add_argument("--max_len", type=float, default=100.0, help="Maximum IBD length")
args = parser.parse_args()

# Save the prefix of the input file  
input_path = Path(args.input_file)
filename = input_path.name.split(".csv")[0]
prefix = filename.split("IBD-Analysis_IBD_")[1]

# Modify the length of IBD segments analyzed
MLE_estim_error.min_len = args.min_len
MLE_estim_error.max_len = args.max_len

### Utility functions ### 

# Euclidean function between two individuals
def pw_dist(x,y):
    """Return pw. distance between two locations"""
    d = np.sqrt(np.sum((x-y)**2))
    return d

# Construct MLE object from IBD + coordinates data
def get_mle_object(df_iids, df_ibds):
	start_list = df_iids[["x", "y"]].values
	l = len(start_list)
	n_pairs = (l * (l - 1) // 2) 
	pair_IBD = [[] for _ in range(n_pairs)]  
	pair_nr = np.ones(n_pairs)
	assert(np.min(pair_nr) >= 1)  
	pair_dist = [pw_dist(start_list[i], start_list[j]) for i in range(l) for j in range(i)]
	id_to_index = {id_: idx for idx, id_ in enumerate(df_iids["Unnamed: 0"])}
	for _, row in df_ibds.iterrows():
	    ibd_length = row["length"]
	    try:
		    ind1 = id_to_index[row["individual1"]]
		    ind2 = id_to_index[row["individual2"]]
	    except KeyError as e:
		    print(f"❌ ID not found in df_iids: {e}")
		    continue
	    j, i = min(ind1, ind2), max(ind1, ind2) 
	    if i != j:							 
		    pair_IBD[(i * (i - 1)) // 2 + j].append(ibd_length)  
	assert(len(pair_dist) == len(pair_IBD))
	assert(len(pair_dist) == len(pair_nr))
	pair_dist, pair_IBD, pair_nr = np.array(pair_dist), np.array(pair_IBD, dtype=object), np.array(pair_nr)
	mle_analyze = MLE_analyse(0, pair_dist, pair_IBD, pair_nr, position_list=start_list, error_model=True)
	return mle_analyze

# Run MLE analysis 
def run_mle_analyze(df_iids, df_ibds, model=args.model, start_params=args.start_params, g=args.g): 		 
    mle_analyze = get_mle_object(df_iids, df_ibds)
    mle_analyze.create_mle_model(model=model, g=g, start_param=np.array(start_params), diploid=0)
    mle_analyze.mle_analysis_error()  
    
    output_txt = os.path.join(args.output_folder.rstrip("/") + '/Results_MLE_' + prefix + '.txt')
    output_csv = os.path.join(args.output_folder.rstrip("/") + '/Results_MLE_' + prefix + '.csv')

    with open(output_txt, "w") as f:
        with redirect_stdout(f):
            mle_analyze.mle_analysis_error()
    
    df_results = parse_mle_output_txt(output_txt, output_csv)
    print("✅ Numerical results extracted from the text output")

# Extract parameters from output text 
def parse_mle_output_txt(output_txt_path, output_csv_path):
    param_rename_map = {"x1": "densité", "const": "dispersion"}
    with open(output_txt_path, "r") as f:
        lines = f.readlines()
    results = []
    start_parsing = False
    for line in lines:
        # Detect start of coefficients table 
        if re.match(r"\s*coef\s+std err\s+z\s+P>\|z\|\s+\[0.025\s+0.975\]", line):
            start_parsing = True
            continue
        # If it's a empty or a decorative ligne, stop 
        if start_parsing and (line.strip() == "" or "=" in line):
            break
        # Read table lines
        if start_parsing:
            parts = line.strip().split()
            if len(parts) == 7:
                param, coef, std_err, z, pval, ci_low, ci_high = parts
                param_rename = param_rename_map.get(param, param)
                results.append({"parameter": param_rename, "coef": float(coef), "std_err":float(std_err), "z": float(z), "pval": float(pval), "lower_CI_95%": float(ci_low), "upper_CI_95%": float(ci_high)})
    # Create DataFrame
    df = pd.DataFrame(results)
    # Save into a .csv 
    df.to_csv(output_csv_path, index=False)
    return df

### Run MLE analysis ### 
df_ibds = pd.read_csv(args.input_file, index_col=0)
df_iids = pd.read_csv(args.loc_file)

run_mle_analyze(df_iids, df_ibds, start_params=args.start_params, g=args.g, model=args.model)
