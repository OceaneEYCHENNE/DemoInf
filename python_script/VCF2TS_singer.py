import os
import subprocess
from pathlib import Path

### Parameters ### 
singer_master = "/mnt/DATA1/VCF/SINGER/singer-0.1.8-beta-linux-x86_64/singer_master"
merge_script = "/mnt/DATA1/VCF/SINGER/singer-0.1.8-beta-linux-x86_64/merge_ARG.py"
vcf_template = "/mnt/DATA1/VCF/VCF_generated_from_TS/forw_homo_m1.0_N1.0_theo_rep{rep}_subsample"
output_base_dir = "/mnt/DATA1/VCF/SINGER/ARG"

mutation_rate = 1e-8
chunk_size = 5_000_000  # length of genomics segments (5Mb par chunk)
chrom_length = 100_000_000  # chromosome total length 
reps = 10  # number of replicates

### Utility functions ### 

# Execute SINGER on one chunk #
def run_chunk(rep, start, end, output_dir):
    output_prefix = output_dir / f"rep{rep}_start{int(start)}_end{int(end)}"
    vcf_file = vcf_template.format(rep=rep)
    cmd = [
        singer_master,
        "-m", str(mutation_rate),
        "-vcf", str(vcf_file),
        "-output", str(output_prefix),
        "-start", str(start),
        "-end", str(end)
    ]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

# Construct file table for merged ARG step ### 
def create_file_table(rep, chunk_files, file_table_path):
    for node_file, branch_file, mut_file, start in chunk_files:
        for f in [node_file, branch_file, mut_file]:
            if not f.exists():
                raise FileNotFoundError(f"No find chunk file {f} for replicate {rep} !")
    with open(file_table_path, "w") as f:
        for i, (node_file, branch_file, mut_file, start) in enumerate(chunk_files):
            f.write(f"{node_file} {branch_file} {mut_file} {start}\n")


### Loop on all replicates ### 
for rep in range(1, reps + 1):
    print(f"\n=== Simulation rep {rep} ===")
    
    # Create an output file per replicate
    rep_dir = Path(output_base_dir) / f"rep{rep}"
    rep_dir.mkdir(parents=True, exist_ok=True)
    
    # Chunk generation
    chunk_files = []
    for start in range(0, chrom_length, chunk_size):
        end = min(start + chunk_size, chrom_length)
        output_prefix = rep_dir / f"rep{rep}_start{start}_end{end}"
        node_file = output_prefix.with_name(output_prefix.name + "_nodes_0.txt")
        branch_file = output_prefix.with_name(output_prefix.name + "_branches_0.txt")
        mut_file = output_prefix.with_name(output_prefix.name + "_muts_0.txt")
        chunk_files.append((node_file, branch_file, mut_file, start))
        if not (node_file.exists() and branch_file.exists() and mut_file.exists()):
            print(f"Chunk {start}-{end} missing, lauching de run_chunk")
            run_chunk(rep, start, end, rep_dir)
        else:
            print(f"Chunk {start}-{end} has already been generated, so we skip run_chunk")
            
    # Create file_table
    file_table_path = rep_dir / "file_table.txt"
    create_file_table(rep, chunk_files, file_table_path)
    
    # Final merge of ARGs
    merged_output = rep_dir / f"rep{rep}_merged_ARG.trees"
    if not merged_output.exists():
        cmd_merge = [
            "python", merge_script,
            "--file_table", str(file_table_path),
            "--output", str(merged_output)
        ]
        print("\nMerging chunks for rep", rep)
        subprocess.run(cmd_merge, check=True)
        print(f"Done rep {rep}, merged ARG saved to {merged_output}\n")
    else:
        print(f"Merge déjà fait pour rep {rep}, fichier {merged_output} existe")
