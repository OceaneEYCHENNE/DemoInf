import subprocess
import sys
import os
import pandas as pd 
import csv
import argparse
import math
import re

### Define arguments ###
# Initialize parser
parser = argparse.ArgumentParser(description="Run plotmaps on MAPS estimation")

# Add arguments with a default value or not 
parser.add_argument("--input_file", "-i", type=str, required=True, help="Path to the folder which contains MAPS Results")
parser.add_argument("--output_folder", "-o", type=str, required=True, help="Path for output log")
args = parser.parse_args()

# Save the prefix of the input file  
rep_folder = os.path.basename(os.path.dirname(args.input_file))
prefix = rep_folder
condition = re.sub(r"_rep\d+$", "", rep_folder)

# Create output folder for each condition
main_output_dir = os.path.join(args.output_folder, "MAPS_RESULTS")
condition_dir = os.path.join(main_output_dir, condition)
os.makedirs(condition_dir, exist_ok=True)


# Function to execute a R script 
def run_r_script(r_script_path, args=None, output_txt=None):
    import subprocess, sys
    try:
        command = ["Rscript", r_script_path]
        if args:
            command.extend(map(str, args))
        print(f"\n>>> Running R script: {' '.join(command)}")
        if output_txt:
            with open(output_txt, "w") as f:
                subprocess.run(command, check=True, stdout=f, stderr=subprocess.STDOUT)
        else:
            subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"❌ Error executing R script: {e}")
        sys.exit(1)

# Extract statistics from log
def extract_stats_from_log(log_path, csv_m_path, csv_N_path):
    with open(log_path, "r") as f:
        lines = [line.strip() for line in f.readlines()]
    
    # Highlight a block of text using a keyword
    def extract_block(keyword):
        for i, line in enumerate(lines):
            if keyword in line:
                header = lines[i+1]
                values = lines[i+2]
                return header, values
        raise ValueError(f"Keyword '{keyword}' not found")

    # Reconstruct correctly column names and maps numerical values
    def parse_line(header_line, values_line):
        raw_headers = header_line.split()
        values = values_line.split()

        headers = []
        i = 0
        while i < len(raw_headers):
            if raw_headers[i] in ["1st", "3rd"]:
                headers.append(raw_headers[i] + " " + raw_headers[i+1])  # "1st Qu."
                i += 2
            else:
                headers.append(raw_headers[i])
                i += 1
        return headers, values

    # Extract migration data
    m_header, m_values = extract_block("plotting migration surface")
    migration_names, migration_values = parse_line(m_header, m_values)

    # Extract population size data
    N_header, N_values = extract_block("plotting population-size surface")
    popsize_names, popsize_values = parse_line(N_header, N_values)

    # Write CSV migration
    with open(csv_m_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["stat", "migration_surface"])
        for name, val in zip(migration_names, migration_values):
            writer.writerow([name, val])

    # Write CSV population size
    with open(csv_N_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["stat", "population_size_surface"])
        for name, val in zip(popsize_names, popsize_values):
            writer.writerow([name, val])

# Extract value of min / max / mean 
def extract_min_max(csv_path):
    df = pd.read_csv(csv_path, sep="\t") 
    min_vals = {}
    max_vals = {}
    mean_vals = {}
    for col in df.columns:
        if col == "stat":
            continue
        min_val = df.loc[df["stat"] == "Min.", col].values[0]
        max_val = df.loc[df["stat"] == "Max.", col].values[0]
        mean_val = df.loc[df["stat"] == "Mean", col].values[0]
        min_vals[col] = float(min_val)
        max_vals[col] = float(max_val)
        mean_vals[col] = float(mean_val)
    return min_vals, max_vals, mean_vals


## Script path ##
plotmaps = "/mnt/DATA1/SCRIPT_R/plotmaps.R"

## Run plotmaps a first time for obtain the summary tables ##
output_results = os.path.join(condition_dir, f"Results_MAPS_{prefix}.txt")
m_min = 0.000001    # nominal value
m_max = 100         # nominal value
N_min = 0.01          # nominal value
N_max = 5000        # nominal value
run_r_script(plotmaps, args=[args.input_file, condition_dir, prefix, m_min, m_max, N_min, N_max], output_txt = output_results)

# Extract the summary stats into a log file separated 
csv_m_path = os.path.join(condition_dir,  f"summary_m_{prefix}.csv")
csv_N_path = os.path.join(condition_dir,  f"summary_N_{prefix}.csv")
extract_stats_from_log(output_results, csv_m_path, csv_N_path)

# Bind migration stats and population size stat into a summary file 
df_m = pd.read_csv(csv_m_path)
df_N = pd.read_csv(csv_N_path)
csv_summary = pd.DataFrame()
csv_summary["stat"] = df_m["stat"]
csv_summary["migration_surface"] = df_m.iloc[:, 1]
csv_summary["population_size_surface"] = df_N.iloc[:, 1]
summary_path = os.path.join(condition_dir, f"Summary_{prefix}.csv")
csv_summary.to_csv(summary_path, index=False, sep="\t")


## Run a new time plotmaps for plot figure with representative scale  ## 
# Extract minminal and maximal value of summary stats 
min_vals, max_vals, mean_vals = extract_min_max(summary_path)
print("Min values:", min_vals)
print("Max values:", max_vals)
print("Mean values:", mean_vals)
m_min = min_vals["migration_surface"]
m_max = max_vals["migration_surface"]
N_min = min_vals["population_size_surface"]
N_max = max_vals["population_size_surface"]

# Round these values for represent better gradient 
m_min = max(math.floor(m_min / 10) * 10, 1e-6)
m_max = math.ceil(m_max / 100) * 100
N_min = max(math.floor(N_min / 10) * 10, 1e-6)
N_max = math.ceil(N_max / 1000) * 1000
print(f"🔁 Re-run avec m ∈ [{m_min}, {m_max}], N ∈ [{N_min}, {N_max}]")

# Re-run plotmaps 
run_r_script(plotmaps, args=[args.input_file, condition_dir, prefix, m_min, m_max, N_min, N_max], output_txt=None)

