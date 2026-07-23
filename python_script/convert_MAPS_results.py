import os
import pandas as pd
import numpy as np
import re
import glob

### Load parameters ###
input_dir = "/mnt/DATA1/MAPS_RESULTS"
output_dir = "/mnt/DATA1/MAPS_RESULTS_KM"

origin_lon = 3.0
origin_lat = 42.5

km_per_deg_lat = 111.0
km_per_deg_lon = 111.0 * np.cos(np.radians(origin_lat))

os.makedirs(output_dir, exist_ok=True)

### Extract true value ### 
def extract_true_value(prefix, param_type):
    if param_type == "sigma":
        match = re.search(r"m([0-9.]+)", prefix)
    elif param_type in ["De", "Ne"]:
        match = re.search(r"N([0-9.]+)", prefix)
    else:
        return np.nan

    return float(match.group(1)) if match else np.nan

### Convert coord lon/lat -> km ###
def lonlat_to_km(df):
    df = df.copy()
    
    x_km = (df["x"] - origin_lon) * km_per_deg_lon
    y_km = (df["y"] - origin_lat) * km_per_deg_lat

    # spatial discretisation
    df["x_km"] = np.floor(x_km).astype(int)
    df["y_km"] = np.ceil(y_km).astype(int)

    return df

### Load MAPS files ### 
files = glob.glob(os.path.join(input_dir, "**", "MAPS_*.csv"), recursive=True)

print(f"Found {len(files)} files")

### Convert coordinates and add theoritical values ### 

for f in files:

    name = os.path.basename(f)
    condition = os.path.basename(os.path.dirname(f))

    df = pd.read_csv(f)

    # identify parameters
    if "sigma" in name:
        param_type = "sigma"
    elif "Ne" in name or "De" in name:
        param_type = "De"
    else:
        continue

    # convert coordinates
    df = lonlat_to_km(df)

    # add true value
    true_val = extract_true_value(condition, param_type)
    df["true_value"] = true_val

    out_subdir = os.path.join(output_dir, condition)
    os.makedirs(out_subdir, exist_ok=True)

    out_path = os.path.join(out_subdir, name.replace(".csv", "_km.csv"))

    df.to_csv(out_path, index=False)

    print(f"✔ {name} | condition={condition} | true={true_val}")

print("\n✅ Done")


