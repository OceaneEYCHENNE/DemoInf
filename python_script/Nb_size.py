import os
import math
import numpy as np
import pandas as pd
import tskit
import pyslim
import re
import argparse
from scipy.spatial import cKDTree
from pathlib import Path 

### Define arguments ###
# Initialize parser
parser = argparse.ArgumentParser(description="Calculate neighbourhood size Nb and effective density De from tree sequence files")

# Add arguments with a default value or not 
parser.add_argument("--input_TS", "-i", type=str, required=True, help="Path to the folder which contains tree-sequences")
parser.add_argument("--output_folder", "-o", type=str, required=True, help="Path for final output csv")
args = parser.parse_args()

# Save the prefix of the input file  
input_path = Path(args.input_TS)
prefix = input_path.name.split(".trees")[0]

### Utility functions ### 

# Extraction a parameter (sigma) from the file name 
def extract_sigma(sim_name):
        match = re.search(r"m([0-9]*\.?[0-9]+)", sim_name)
        if match:
            return float(match.group(1))
        else:
            raise ValueError(f"Unable to extract “sigma” from the name '{sim_name}'")

# Calculate the average number of neighbours within a given spatial radius 
def mean_neighbors_in_radius(ts_path,
                             sigma = None,
                             alive_only=True,
                             generation=0, 
                             include_self=False):
    # Define spatial radius (linked to sigma)
    if sigma is None:
        sigma = extract_sigma(ts_path.split("/")[-1]) # prend juste le nom du fichier 
    radius = sqrt(2.0) * sigma
    
    # Load tree sequence
    ts = tskit.load(ts_path)
    
    # Select individuals (alive or all)
    if alive_only:
        inds = pyslim.individuals_alive_at(ts, generation)
        if len(inds) == 0:
            raise ValueError(f"No living individuals found in that generation {generation}")
    else:
        inds = np.arange(ts.num_individuals, dtype=int)
        
    # Extract spatial coordinates
    locs = ts.individuals_location[inds, :2] 
    
    # Search for neighbours using a spatial tree (KD-tree)
    tree = cKDTree(locs)
    neighbors = tree.query_ball_point(locs, r=radius)
    
    # Count neighbours
    counts = np.array([len(lst) for lst in neighbors], dtype=int)
    if not include_self:
        counts -= 1
    counts[counts < 0] = 0
    mean_neighbors = counts.mean()
    return mean_neighbors

# Convert Nb into effective density (De)
def add_De_from_Nb(df, sim_col='sim', Nb_col='Nb'):
    df = df.copy()
    df["sigma"] = df[sim_col].apply(extract_sigma)
    df["De"] = df[Nb_col] / (4 * math.pi * df["sigma"]**2)
    return df

### Read all trees sequences ###
folder_path = args.input_TS
tree_files = [f for f in os.listdir(folder_path) if f.endswith(".trees")]

### Calculate spatial statistics ### 
results = []

for f in tree_files:
    ts_path = os.path.join(folder_path, f)
    try:
        nb = mean_neighbors_in_radius(ts_path)
        results.append({"sim": f, "Nb": nb})
    except Exception as e:
        print(f"Erreur pour {f} : {e}")
        results.append({"sim": f, "Nb": np.nan})

df = pd.DataFrame(results)
print(df)

# Calculate De 
df_De = add_De_from_Nb(df)

# Save results into a csv file
df.to_csv(os.path.join(args.output_folder.rstrip("/") + "/mean_neighbors_summary.csv"), index=False)
df_De.to_csv(os.path.join(args.output_folder.rstrip("/") + "/mean_neighbors_De_summary.csv"), index=False)

