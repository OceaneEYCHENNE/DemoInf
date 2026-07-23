import argparse
import tskit 
import tsinfer
import tsdate
import cyvcf2
import numpy as np
import pandas as pd 
import re
import csv
from pathlib import Path
import json

### Define arguments ###
# Initialize parser
parser = argparse.ArgumentParser(description="Script to convert VCF file into a tree sequence with tsinfer/tsdate algorithm")

# Add arguments with a default value or not 
parser.add_argument("--input_file", "-i", type=str,required=True, help="Path for input vcf")
parser.add_argument("--output_folder", "-o", type=str, required=True, help="Path for output")
args = parser.parse_args()

# Save the prefix of the input file  
input_path = Path(args.input_file)
prefix = input_path.name.split(".vcf")[0]

### Construct sample data ###
vcf_samples = cyvcf2.VCF(args.input_file).samples
    
samples_file = args.output_folder.rstrip("/") + "/" + prefix + "_tsinfer.samples"
L = 100000000 

with tsinfer.SampleData(path=str(samples_file), sequence_length=L) as samples:
    # Add diploid individuals
    vcf_ids = vcf_samples
    for vcf_name in vcf_ids:
        uid = f"{vcf_name}"
        samples.add_individual(ploidy=2, metadata={"uid": uid})
    # Add variants (SNPs only) 
    seen_positions = set() 
    for variant in cyvcf2.VCF(args.input_file):
        if not variant.is_snp:
            continue
        # management of duplicate positions
        pos = variant.POS
        epsilon = 1e-6
        while pos in seen_positions:
            pos += epsilon 
        seen_positions.add(pos)
        # extract genotypes and alleles
        gts = variant.genotype.array()[:, :2].flatten()
        # add sites 
        samples.add_site(position=pos, genotypes=gts, alleles=[variant.REF] + list(variant.ALT))

### Construct mapping file ### 
mapping_file = args.output_folder.rstrip("/") + "/" + prefix + "_ind_map.csv"

with open(mapping_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["tsinfer_id", "vcf_name"])
    
    for ind in samples.individuals():
        writer.writerow([ind.id, ind.metadata["uid"]])

## Inference of the tree sequence (tsinfer)
ts_inferred = tsinfer.infer(samples, num_threads=4, progress_monitor=True)
ts_inferred.dump(args.output_folder.rstrip("/") + "/" + prefix + "_inferred.trees")

## Datation of the tree sequence (tsdate)
ts_simplify = ts_inferred.simplify(keep_unary = False)
ts_dated = tsdate.date(ts_simplify, mutation_rate=1e-8)
ts_dated.dump(args.output_folder.rstrip("/") + "/" + prefix + "_dated.trees")

## Keep metadata individual information 
tables = ts_dated.dump_tables()
for i, ind in enumerate(tables.individuals):
    tables.individuals[i] = ind.replace(
        metadata=json.dumps({"name": vcf_samples[i]}).encode()
    )
ts_dated_with_meta = tables.tree_sequence()
ts_dated_with_meta.dump(args.output_folder.rstrip("/") + "/" + prefix + "_dated_metadata.trees")
