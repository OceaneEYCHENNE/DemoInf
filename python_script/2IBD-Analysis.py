import argparse  
import numpy as np
import pandas as pd 
from pathlib import Path

### Define arguments ###
# Initialize parser
parser = argparse.ArgumentParser(description="Script to convert IBD segments file and coordinates file into the valid format for demographic inference by IBD-Analysis")

# Add arguments with a default value or not 
parser.add_argument("--IBD_file", "-i", type=str,required=True, help="Path for IBD file")
parser.add_argument("--output_folder", "-o", type=str, required=True, help="Path for output ")
parser.add_argument("--locs_file", "-loc",type=str, help="Path for location file")
args = parser.parse_args()
 
# Save the prefix of the input file  
input_path = Path(args.IBD_file)
filename = input_path.name.split(".csv")[0]
prefix = filename.split("IBD_coord_")[1]
 
### Import files ###
dfIBD = pd.read_csv(args.IBD_file)
print(f"Loaded a csv file containing {dfIBD.Pair.count()} IBD segments")

dfLOC = pd.read_csv(args.locs_file).values
print(f"Loaded a csv file containing {len(dfLOC)} individuals locations")
 
### Convert IBD file format ###
 
# Keep only the columns of individual number, x-coordinate and y-coordinate
dfIBD_filtered = dfIBD[['Ind1', 'Ind2', 'Left', 'Right', 'Length']]

# Add a new column with the number of IBD segments 
dfIBD_filtered.insert(0, '', range(len(dfIBD_filtered))) 

# Rename column header
dfIBD_filtered.columns= ['', 'individual1', 'individual2', 'start', 'endpoint', 'length']

# Convert column start, endpoint and length into cM 
dfIBD_filtered.loc[:, 'length'] = dfIBD_filtered['length'].astype(float) * 1e-6
dfIBD_filtered.loc[:, 'start'] = dfIBD_filtered['start'].astype(float) * 1e-6
dfIBD_filtered.loc[:, 'endpoint'] = dfIBD_filtered['endpoint'].astype(float) * 1e-6

# Export new IBD file 
dfIBD_filtered.to_csv(args.output_folder.rstrip("/") + '/IBD-Analysis_IBD_' + prefix + '.csv', index=False)
print(f"Converted IBD file into the valid format for IBD-Analysis")


### Convert locations file format ###
# Keep only the columns of individual number, x-coordinate and y-coordinate
dfLOC_filtered = dfLOC[:, [0, 2, 3]].copy()

# Remove column header 0
dfLOC_filtered = pd.DataFrame(dfLOC_filtered, columns=['', 'x', 'y'])
dfLOC_filtered[''] = dfLOC_filtered[''].astype(int)

# Export new locations file 
dfLOC_filtered.to_csv(args.output_folder.rstrip("/") + '/IBD-Analysis_Coord_' + prefix + '.csv', index=False)
print(f"Converted locations file into the valid format for IBD-Analysis")
 
