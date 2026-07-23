import argparse
import msprime
import random
import pandas as pd
import numpy as np
import tskit
import gzip
from pathlib import Path

### Define arguments ###
# Initialize parser
parser = argparse.ArgumentParser(description="Convert tree sequences files into VCF files ")

# Add arguments with a default value or not 
parser.add_argument("--input_file", "-i", type=str,required=True, help="Path for input ts file")
parser.add_argument("--output_folder", "-o", type=str, required=True, help="Path for output folder")
args = parser.parse_args()

# Save the prefix of the input file  
input_path = Path(args.input_file)
filename = input_path.name.split(".trees")[0]
prefix = filename.split("forw_")[1]
output_folder = Path(args.output_folder)

### Import tree-sequence of an SLiM simulation ###
ts = tskit.load(args.input_file)
print(f"Loaded a simulated tree sequence of {ts.num_samples} haploid genomes from {ts.num_individuals} diploid individuals:")
print(f"containing {ts.num_sites} variable sites over {ts.sequence_length} base pairs")

### Subsample and simplify the simulated tree sequence for keep only alive individuals
alive = [ind.id for ind in ts.individuals() if ind.time == 0]
samples = []
for ind in alive:
    samples.extend(ts.individual(ind).nodes)
ts_simplified = ts.simplify(samples=samples, keep_input_roots=False)

### Add mutations on tree sequence ### 
ts_mutated = msprime.sim_mutations(ts_simplified, rate=1e-8, model=msprime.JC69(), random_seed=42)

### Identify individuals with samples ### 
sample_set = set(ts_mutated.samples())
individual_ids = [
    ind.id for ind in ts_mutated.individuals()
    if any(n in sample_set for n in ind.nodes)
]
individual_names = [str(ind_id) for ind_id in individual_ids]

### Save VCF file with SLiM individuals names ###
vcf_path = output_folder / f"{prefix}.vcf"
with open(vcf_path, "w") as f:
    ts_mutated.write_vcf(f, individuals=individual_ids, individual_names=individual_names)
    print(f"VCF written: {vcf_path}")
    
# Write coordinates
rows = [
    {"id": ind_id,
     "x": ts_mutated.individual(ind_id).location[0],
     "y": ts_mutated.individual(ind_id).location[1]}
    for ind_id in individual_ids
]
coords_df = pd.DataFrame(rows)
coords_path = output_folder / f"Coord_{prefix}.csv"
coords_df.to_csv(coords_path, index=False)
print(f"Coordinates written: {coords_path}")

