import pandas as pd
import numpy as np
from pathlib import Path
import argparse  

### Define arguments ###
# Initialize parser
parser = argparse.ArgumentParser(description="Script to convert IBD segments file and coordinates file into the valid format for demographic inference by MAPS")

# Add arguments with a default value or not 
parser.add_argument("--IBD_file", "-i", type=str,required=True, help="Path for IBD file")
parser.add_argument("--output_folder", "-o", type=str, required=True, help="Path for output ")
parser.add_argument("--locs_file", "-loc",type=str, help="Path for location file")
parser.add_argument("--habitat_size", "-W",type=str, default = '100', help="Value of size habitat i.e. square habitat")
parser.add_argument("--lower_bound", "-lwr_b", type=str, default = '2', help="Value of lower bound (in cM)")
parser.add_argument("--upper_bound", "-upp_b", type=str, default = 'Inf', help="Value of upper bound (in cM)")
args = parser.parse_args()

 
# Save the prefix of the input file  
input_path = Path(args.IBD_file)
filename = input_path.name.split(".csv")[0]
prefix = filename.split("IBD_coord_")[1]

### Import files ###
dfIBD = pd.read_csv(args.IBD_file)
print(f"Loaded a csv file containing {dfIBD.Pair.count()} IBD segments")

dfLOC = pd.read_csv(args.locs_file)
print(f"Loaded a csv file containing {len(dfLOC)} individuals locations")
dfLOC = dfLOC.iloc[:, 2:4]  # Keep only columns of coordinates x and coordinates y
dfLOC = pd.DataFrame(np.repeat(dfLOC.values, 2, axis=0), columns=dfLOC.columns)             # Double each line to obtain coordinate values for the 2 haplotypes of each individual


## Prepare IDs ##
ids = list(range(len(dfLOC)))                                                               # Create list of haplotype IDs
n = len(ids)                                                                                # Calculate total number of haplotypes (2 * individuals number)


### Prepare coordinates file ###
# Conversion function 
def convert_km_to_lonlat(x_km, y_km, origin_lon=3.0, origin_lat=42.5):
    """
    Convert km coordinates to geographic coordinates (longitude, latitude).
    """
    km_per_deg_lat = 111.0
    km_per_deg_lon = 111.0 * np.cos(np.radians(origin_lat))

    lon = origin_lon + np.array(x_km) / km_per_deg_lon
    lat = origin_lat + np.array(y_km) / km_per_deg_lat

    return pd.DataFrame({'lon': lon, 'lat': lat})
    
# Convert to latitude/longitude
dfLOC_geo = convert_km_to_lonlat(dfLOC.iloc[:, 0], dfLOC.iloc[:, 1])

# Export new locations file 
np.savetxt(args.output_folder.rstrip("/") + '/' + prefix + '_' + args.lower_bound + '_' + args.upper_bound + '.coord', dfLOC_geo.to_numpy(), fmt='%.6f', delimiter=' ')  
print(f"Converted locations file into the valid format for MAPS")

### Create outer file ### 
outer = [
    (0, 0),
    (0, float(args.habitat_size)),
    (float(args.habitat_size), float(args.habitat_size)),
    (float(args.habitat_size), 0),
    (0, 0)
]
outer = pd.DataFrame(outer)

# Convert to latitude/longitude
outer_geo = convert_km_to_lonlat(outer.iloc[:, 0], outer.iloc[:, 1])
outer_geo = np.array(outer_geo, dtype=float)

output_path = args.output_folder.rstrip("/") + '/' + prefix + '_' + args.lower_bound + '_' + args.upper_bound + '.outer'
with open(output_path, "w") as f:
    for x, y in outer_geo:
        f.write(f"{x:.6f} {y:.6f}\n")
print(f"Saved outer file")


### Create IBD file ### 
## Initiate similarity matrix ## 
ibd_summary = np.zeros((n, n), dtype=int)                                                   # Create an empty n*n matrix in order to store the number of shared IBD segments

## Filter IBD segment by length
lengths = dfIBD['Length'] / 1e6                                                             # Calculate length of each IBD segment and convert into cM   
selected_inds = dfIBD[(lengths > int(args.lower_bound)) & (lengths < float(args.upper_bound))].index      # Select only IBD segments which length is superior to 2 cM

## Filling of similarity matrix  
for i in selected_inds:                                                                     # For each filtered segment
    id1 = dfIBD.at[i, 'Hap1']                                                               # Recover haplotype1 ID 
    id2 = dfIBD.at[i, 'Hap2']                                                               # Recover haplotype2 ID
    if id1 in ids and id2 in ids:                                                           # If this haplotype ID is in the ID list
        ibd_summary[id1, id2] += 1                                                          # Increment the matrix by the number of shared IBD segments of this haplotype pair
        ibd_summary[id2, id1] = ibd_summary[id1, id2]                                       # As relation is symmetric, fill the 2 entries

## Save files 
np.savetxt(args.output_folder.rstrip("/") + '/' + prefix + '_' + args.lower_bound + '_' + args.upper_bound + '.sims', ibd_summary, fmt='%d', delimiter=' ')                 # Save similarity matrix file
print(f"Converted IBD file into the valid format for MAPS")
