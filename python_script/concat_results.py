import pandas as pd
import os
import glob
import argparse
import re
import numpy as np

### Define arguments ###
# Initialize parser
parser = argparse.ArgumentParser(description="Concatenate MAPS and IBD Analussi results (even if MAPS not available")

# Add arguments with a default value or not 
parser.add_argument("--input_MAPS", "-M", type=str, required=True, help="Path to MAPS Results folder")
parser.add_argument("--input_IBD_Analysis", "-I", type=str, required=True, help="Path to IBD-Analysis Results folder")
parser.add_argument("--output_folder", "-o", type=str, required=True, help="Output folder")
args = parser.parse_args()

# Define Folders containing files 
maps_dir = args.input_MAPS
ibd_dir = args.input_IBD_Analysis

### List files ### 
maps_files = {re.search(r"Summary_(.*)\.csv", f).group(1): f
              for f in glob.glob(os.path.join(maps_dir, "Summary_*.csv"))}

ibd_files = {re.search(r"Results_MLE_(.*)\.csv", f).group(1): f
             for f in glob.glob(os.path.join(ibd_dir, "Results_MLE_*.csv"))}

print(f"📁 MAPS files detected: {len(maps_files)}")
print(f"📁 IBD files detected: {len(ibd_files)}")

### Type detection function ### 
def get_sim_type(prefix):
    # Case 1 : prefix start with 'homo'
    if prefix.startswith("homo"):
        # Check simulation types
        lc_match = re.search(r"life_cycle_[^_]+", prefix)
        h_match = re.search(r"herm_[^_]+", prefix)
        l_match = re.search(r"lap", prefix)
        s_match = re.search(r"stu", prefix) 
        temp_m_match = re.search(r"temporal_m", prefix)
        temp_N_match = re.search(r"temporal_N", prefix)
        mate1_match = re.search(r"mate1.5", prefix) 
        mate0_match = re.search(r"mate0.5", prefix) 
        comp1_match = re.search(r"comp1.5", prefix) 
        comp0_match = re.search(r"comp0.5", prefix)
        edge_match = re.search(r"edge", prefix)
        disp_match = re.search(r"disp", prefix)
        sigma_match = re.search(r"sigma", prefix) 
        singer_ARG_match = re.search(r"singer_ARG", prefix)
        singer_match = re.search(r"singer", prefix)
        relate_ARG_match = re.search(r"relate_ARG", prefix)
        ts_ARG_match = re.search(r"tsinfer_ARG", prefix)
        ts_match = re.search(r"tsinfer", prefix)
        Em3000_match = re.search(r"herm_Emarg_3000", prefix)
        Em4000_match = re.search(r"herm_Emarg_4000", prefix)
        Em5000_match = re.search(r"herm_Emarg_5000", prefix)
        Ds300_match = re.search(r"herm_Dsarg_300", prefix)
        Ds400_match = re.search(r"herm_Dsarg_400", prefix)
        Ds500_match = re.search(r"herm_Dsarg_500", prefix)
        
        if lc_match: return lc_match.group(0)
        if h_match: return h_match.group(0)
        if l_match: return l_match.group(0)
        if s_match: return s_match.group(0) 
        if temp_m_match: return temp_m_match.group(0)
        if temp_N_match: return temp_N_match.group(0)
        if mate1_match: return mate1_match.group(0)
        if mate0_match: return mate0_match.group(0)
        if comp1_match: return comp1_match.group(0)
        if comp0_match: return comp0_match.group(0)
        if edge_match: return edge_match.group(0)
        if disp_match: return disp_match.group(0)
        if sigma_match: return sigma_match.group(0)
        if singer_ARG_match: return singer_ARG_match.group(0)
        if singer_match: return singer_match.group(0)
        if relate_ARG_match: return relate_ARG_match.group(0) 
        if ts_ARG_match: return ts_ARG_match.group(0) 
        if ts_match: return ts_match.group(0) 
        if Em3000_match: return Em3000_match.group(0)
        if Em4000_match: return Em4000_match.group(0)
        if Em5000_match: return Em5000_match.group(0)
        if Ds300_match: return Ds300_match.group(0)
        if Ds400_match: return Ds400_match.group(0)
        if Ds500_match: return Ds500_match.group(0)
        
        # Search values in the prefix  
        m_match = re.search(r"m([0-9.]+)", prefix)
        N_match = re.search(r"N([0-9.]+)", prefix)
        K_match = re.search(r"K([0-9]+)", prefix)
        t_match = re.search(r"t([0-9]+)", prefix)
        m_val = float(m_match.group(1)) if m_match else None
        N_val = float(N_match.group(1)) if N_match else None
        K_val = int(K_match.group(1)) if K_match else None
        t_val = int(t_match.group(1)) if t_match else None
        if m_val == 1.0 and N_val == 1.0 and K_val == 200 and t_val == 201:
            return "condi_moyenne"
        if m_match and m_val != 1.0: return f"m{m_val}"
        if N_match and N_val != 1.0: return f"N{N_val}"
        if K_match and K_val != 200: return f"K{K_val}"
        if t_match and t_val != 201: return f"t{t_val}"
        return "homo_unknown"
        
    # Case 2 : prefix start with 'hetero'
    if prefix.startswith("hetero"):
        if prefix.startswith("hetero_m_v2"): return "hetero_m_v2"
        if prefix.startswith("hetero_N_v2"): return "hetero_N_v2"
        if prefix.startswith("hetero_m"): return "hetero_m"
        if prefix.startswith("hetero_N"): return "hetero_N"
        return "hetero"

    # Case 3 : prefix start with 'grad'
    if prefix.startswith("grad"):
        if prefix.startswith("grad_m_v2"): return "grad_m_v2"
        if prefix.startswith("grad_N_v2"): return "grad_N_v2"
        if prefix.startswith("grad_m"): return "grad_m"
        if prefix.startswith("grad_N"): return "grad_N"
        return "grad"

    else:
        return "unknown"


### Merge loop based on IBD files ### 
all_results_detailled = []

for prefix, ibd_path in ibd_files.items():
    sim_type = get_sim_type(prefix)
    maps_path = maps_files.get(prefix, None)

    # IBD Analysis
    try:
        ibd_df = pd.read_csv(ibd_path)
    except Exception as e:
        print(f"❌ Error reading IBD {prefix}: {e}")
        continue

    # Extract IBD parameters
    try:
        disp_row = ibd_df[ibd_df["parameter"] == "dispersion"].iloc[0]
        dens_row = ibd_df[ibd_df["parameter"] == "densité"].iloc[0]
    except Exception as e:
        print(f"❌ Missing IBD rows for {prefix}: {e}")
        continue

    # If MAPS MISSING → MAPS values = NaN
    if maps_path is None:
        print(f"⚠️ MAPS file missing for {prefix}, filling with NaN")
        m_MAPS_stats = {k: np.nan for k in ["Mean", "Min.", "1st Qu.", "3rd Qu.", "Max."]}
        N_MAPS_stats = {k: np.nan for k in ["Mean", "Min.", "1st Qu.", "3rd Qu.", "Max."]}
    else:
        try:
            maps_df = pd.read_csv(maps_path, sep="\t")
            stats = ["Mean", "Min.", "1st Qu.", "3rd Qu.", "Max."]
            maps_stats = maps_df[maps_df["stat"].isin(stats)].set_index("stat")
            m_MAPS_stats = {s: maps_stats.loc[s, "migration_surface"] for s in stats}
            N_MAPS_stats = {s: maps_stats.loc[s, "population_size_surface"] for s in stats}
        except Exception as e:
            print(f"❌ Error reading MAPS {prefix}, filling NaN: {e}")
            m_MAPS_stats = {k: np.nan for k in ["Mean", "Min.", "1st Qu.", "3rd Qu.", "Max."]}
            N_MAPS_stats = {k: np.nan for k in ["Mean", "Min.", "1st Qu.", "3rd Qu.", "Max."]}

    # Append results for m
    all_results_detailled.append({
        "parameter": "m",
        "MAPS_Mean": m_MAPS_stats["Mean"], "MAPS_Min": m_MAPS_stats["Min."],
        "MAPS_Q1": m_MAPS_stats["1st Qu."], "MAPS_Q3": m_MAPS_stats["3rd Qu."],
        "MAPS_Max": m_MAPS_stats["Max."],
        "IBD_coef": disp_row["coef"], "IBD_std_err": disp_row["std_err"],
        "IBD_pval": disp_row["pval"], "IBD_lower95": disp_row["lower_CI_95%"],
        "IBD_upper95": disp_row["upper_CI_95%"],
        "sim": prefix, "type": sim_type
    })

    # Append results for N
    all_results_detailled.append({
        "parameter": "N",
        "MAPS_Mean": N_MAPS_stats["Mean"], "MAPS_Min": N_MAPS_stats["Min."],
        "MAPS_Q1": N_MAPS_stats["1st Qu."], "MAPS_Q3": N_MAPS_stats["3rd Qu."],
        "MAPS_Max": N_MAPS_stats["Max."],
        "IBD_coef": dens_row["coef"], "IBD_std_err": dens_row["std_err"],
        "IBD_pval": dens_row["pval"], "IBD_lower95": dens_row["lower_CI_95%"],
        "IBD_upper95": dens_row["upper_CI_95%"],
        "sim": prefix, "type": sim_type
    })


### Save detailed output ### 
all_results_detailled_df = pd.DataFrame(all_results_detailled)

for sim_type, df_group in all_results_detailled_df.groupby("type"):
    out = os.path.join(args.output_folder, f"Summary_stats_{sim_type}.csv")
    df_group.to_csv(out, index=False)
    print(f"✅ Saved {out}")


### Mean sumary (MAPS may be missing) ### 
all_results_mean = []

for prefix, ibd_path in ibd_files.items():

    sim_type = get_sim_type(prefix)
    maps_path = maps_files.get(prefix, None)

    try:
        ibd_df = pd.read_csv(ibd_path)
    except:
        continue

    # Extract IBD coef and 95% CI 
    disp_row = ibd_df.loc[ibd_df["parameter"] == "dispersion"].iloc[0]
    dens_row = ibd_df.loc[ibd_df["parameter"] == "densité"].iloc[0]
    
    m_IBD = disp_row["coef"]
    m_lower = disp_row["lower_CI_95%"]
    m_upper = disp_row["upper_CI_95%"]

    N_IBD = dens_row["coef"]
    N_lower = dens_row["lower_CI_95%"]
    N_upper = dens_row["upper_CI_95%"]

    # Try reading MAPS mean or fill NaN
    stats = ["Mean", "Min.", "1st Qu.", "3rd Qu.", "Max."]

    if maps_path and os.path.exists(maps_path):
        try:
            maps_df = pd.read_csv(maps_path, sep="\t")
            maps_stats = maps_df[maps_df["stat"].isin(stats)].set_index("stat")

            m_MAPS_stats = {s: maps_stats.loc[s, "migration_surface"] for s in stats}
            N_MAPS_stats = {s: maps_stats.loc[s, "population_size_surface"] for s in stats}

        except:
            m_MAPS_stats = {k: np.nan for k in stats}
            N_MAPS_stats = {k: np.nan for k in stats}
    else:
        m_MAPS_stats = {k: np.nan for k in stats}
        N_MAPS_stats = {k: np.nan for k in stats}


    all_results_mean.append({
    "parameter": "m",
    "MAPS_Mean": m_MAPS_stats["Mean"],
    "MAPS_Min": m_MAPS_stats["Min."],
    "MAPS_Q1": m_MAPS_stats["1st Qu."],
    "MAPS_Q3": m_MAPS_stats["3rd Qu."],
    "MAPS_Max": m_MAPS_stats["Max."],
    "IBD-Analysis": m_IBD,
    "IBD-lower95": m_lower,
    "IBD-upper95": m_upper,
    "sim": prefix,
    "type": sim_type
})

    all_results_mean.append({
    "parameter": "N",
    "MAPS_Mean": N_MAPS_stats["Mean"],
    "MAPS_Min": N_MAPS_stats["Min."],
    "MAPS_Q1": N_MAPS_stats["1st Qu."],
    "MAPS_Q3": N_MAPS_stats["3rd Qu."],
    "MAPS_Max": N_MAPS_stats["Max."],
    "IBD-Analysis": N_IBD,
    "IBD-lower95": N_lower,
    "IBD-upper95": N_upper,
    "sim": prefix,
    "type": sim_type
})

# Save final summary
out = os.path.join(args.output_folder, "Summary_stats.csv")
pd.DataFrame(all_results_mean).to_csv(out, index=False)
print(f"✅ Saved {out}")
