import pyslim
import os
import tskit
import numpy as np
import argparse
from pathlib import Path
import pandas as pd

### Define arguments ###
# Initialize parser
parser = argparse.ArgumentParser(description="Compute effective disperdal sigma_eff from tree sequence files")

# Add arguments with a default value or not 
parser.add_argument("--input_folder", "-i", type=str,required=True, help="Path for input tree-sequence folder")
parser.add_argument("--output_folder", "-o", type=str, required=True, help="Path for output ")
args = parser.parse_args()

### Effective dispersal calculation function ###  

def compute_sigma_eff(ts, sample_only_alive=True):
    displacements = []
    individuals = ts.individuals()
    
    # Filter of living individuals
    if sample_only_alive:
        alive_ids = set(pyslim.individuals_alive_at(ts, 0))
        
    # Loop on individuals and their parents
    for ind in individuals:
        if ind.location is None or len(ind.location) < 2:
            continue
        if sample_only_alive and ind.id not in alive_ids:
            continue
        child_loc = np.array(ind.location[:2])
        for parent_id in ind.parents:
            if parent_id == tskit.NULL:
                continue
            parent = ts.individual(parent_id)
            if parent.location is None or len(parent.location) < 2:
                continue
            parent_loc = np.array(parent.location[:2])
            # quadratic distance parent-childrens
            delta = child_loc - parent_loc
            sq_disp = np.sum(delta ** 2)
            displacements.append(sq_disp)
    displacements = np.array(displacements)
    if len(displacements) == 0:
        raise ValueError("No valid parent-offspring dispersal events found.")
    
    # Estimation of sigma_eff
    sigma_eff = np.sqrt(0.5 * np.mean(displacements))
    return sigma_eff, displacements

### Read all trees sequences ###
folder_path = Path(args.input_folder)
tree_files = [f for f in os.listdir(folder_path) if f.endswith(".trees")]

### Calculate spatial statistics ### 
results = []

for f in tree_files:
    ts_path = os.path.join(folder_path, f)
    try:
        ts = tskit.load(ts_path)
        sigma_eff, sq_disp = compute_sigma_eff(ts, sample_only_alive=True)
        results.append({"sim": f, "sigma_eff": sigma_eff})
    except Exception as e:
        print(f"Erreur pour {f} : {e}")
        results.append({"sim": f, "sigma_eff": np.nan})

df = pd.DataFrame(results)

# Save results into a csv file
df.to_csv(args.output_folder.rstrip("/") + '/sigma_eff.csv', index=False)
