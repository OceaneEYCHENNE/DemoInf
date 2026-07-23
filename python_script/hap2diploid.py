import tskit
import gzip
import pandas as pd 
import numpy as np

### Read individuals in VCF files ### 
vcf_path = "/mnt/DATA1/VCF/VCF_generated_from_ARG_subsample/singletons/forw_homo_m1.0_N1.0_theo_rep1_filtered_singletons.vcf"

opener = gzip.open if vcf_path.endswith(".gz") else open
with opener(vcf_path, "rt") as f:
    for line in f:
        if line.startswith("#CHROM"):
            vcf_samples = line.strip().split("\t")[9:]
            break
print(f"Individuals in the VCF : {len(vcf_samples)}")  


### Load tree sequences ###
ts = tskit.load("/mnt/DATA1/VCF/VCF_generated_from_ARG_subsample/singletons/forw_homo_m1.0_N1.0_theo_rep1.trees")
print(f"Samples dans Relate : {ts.num_samples}")       

assert ts.num_samples == len(vcf_samples) * 2, \
    "The number of Relate samples must be double the number of individuals in the VCF file !"


### Read spatial coordinates ### 
coords_df = pd.read_csv(
    "/mnt/DATA1/VCF/VCF_generated_from_ARG_subsample/Coord_homo_m1.0_N1.0_theo_rep1.csv",
    sep=","
)
coords_df["id"] = coords_df["id"].astype(str)
coords_dict = coords_df.set_index("id")[["x", "y"]].to_dict("index")

### Reconstruct diploid individuals ### 
tables = ts.dump_tables()
tables.individuals.clear()

samples = ts.samples() 
nodes_individual = np.full(tables.nodes.num_rows, -1, dtype=np.int32)

for i, name in enumerate(vcf_samples):
    # Get coordinates
    coords = coords_dict.get(name, {"x": float("nan"), "y": float("nan")})
    
    # Create tskit individuals 
    ind_id = tables.individuals.add_row(location=[coords["x"], coords["y"]], metadata=name.encode())
    
    # Associate haplotypes to individuals
    node_hap1 = samples[2 * i]
    node_hap2 = samples[2 * i + 1]
    nodes_individual[node_hap1] = ind_id
    nodes_individual[node_hap2] = ind_id

tables.nodes.individual = nodes_individual
ts_final = tables.tree_sequence()

### Check and save files ### 
print(f"Individus : {ts_final.num_individuals}")
for ind in list(ts_final.individuals())[:3]:
    name = ind.metadata.decode()
    found = name in coords_dict
    print(f"  {name} → location {ind.location} → nodes {ind.nodes} (coords trouvées: {found})")

ts_final.dump("/mnt/DATA1/VCF/VCF_generated_from_ARG_subsample/singletons/forw_homo_m1.0_N1.0_Relate_ARG_diploid_coord_rep1.trees")
print("Saved.")
