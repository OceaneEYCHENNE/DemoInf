import argparse  
import tskit
import pyslim
import random
import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
import datetime 
import math
import re
import seaborn as sns
from pathlib import Path

### Define arguments ###
# Initialize parser
parser = argparse.ArgumentParser(description="Script to calculate the proportion of genome shared in IBD, as well as the geographical distance to test isolation by distance")

# Add arguments with a default value or not 
parser.add_argument("--input_file", "-i", type=str,required=True, help="Path for IBD file with coordinates informations")
parser.add_argument("--output_folder", "-o", type=str, required=True, help="Path for output folder")
parser.add_argument("--chr_length", "-l", type=int, default = 1e8, help="Length of simulated chromosome")
parser.add_argument("--loc_file", "-loc", type=str,required=True, help="Path for individuals geographics coordinates file")
parser.add_argument("--distrib", "-d", action="store_true", default = True, help="Plot the proportion of genome in shared IBD segments as a function of geographic distance")
args = parser.parse_args()

# Save the prefix of the input file  
input_path = Path(args.input_file)
filename = input_path.name.split(".csv")[0]
prefix = filename.split("IBD_coord_")[1]


### Import files ###
dfIBD = pd.read_csv(args.input_file)
print(f"Loaded a csv file containing {dfIBD.Pair.count()} IBD segments")
dfLOC = pd.read_csv(args.loc_file)
print(f"Recorded {len(dfLOC.Ind)} individual locations")


### Retrieve variable values ###
L = args.chr_length
max_time = int(dfIBD["Time"].max())


### Cut Time by Interval of 5 generations ###
dfIBD['Interval'] = pd.cut(dfIBD['Time'], bins=range(0, max_time, 5), right = False)


### Calculate Geographic Distance between individuals in pair ### 
dfIBD['Dist_Geo'] = np.sqrt((dfIBD['X2'] - dfIBD['X1'])**2 + (dfIBD['Y2'] - dfIBD['Y1'])**2)


### Calculate Genome Proportion shared by IBD segments ###

# Calculate length cumulated for each Pair and each Interval 
Cumulated_length = dfIBD.groupby(['Pair', 'Interval'], observed=False)['Length'].sum().reset_index()

# Replace NaN value by 0
Cumulated_length['Length'] = Cumulated_length['Length'].fillna(0.0)

# Calculate Genome proportion 
Cumulated_length['Prop'] = Cumulated_length['Length'] / L

# Calculate the sum of length and proportion for each pair 
Length_Prop = Cumulated_length.groupby('Pair', observed=False).agg(Length=('Length', 'sum'), Prop=('Prop', 'sum')).reset_index()

# Calculate the mean geographic distance of each pair 
Dist_Geo = dfIBD.groupby('Pair', observed=False)['Dist_Geo'].mean().reset_index()

# Merge the geographic distance with the proprotion of genome shared by Pair 
Length_Prop = Length_Prop.merge(Dist_Geo, on='Pair', how='left')

# Export dataframe into a .csv file #
Length_Prop.to_csv(args.output_folder.rstrip("/") + '/Length_Prop_' + prefix + '.csv', index=False)


### Make the representation of proportion of genome in shared IBD segments as a function of geographic distance ### 
if args.distrib:		
	plt.figure(figsize=(8,6))
	sns.scatterplot(data=Length_Prop, x="Dist_Geo", y="Prop", color="blue", alpha=0.6)
	sns.regplot(data=Length_Prop, x="Dist_Geo", y="Prop", scatter=False, color="red", ci=None)  # Add linear regression
	plt.ylim(-0.05, 1.05)
	plt.xlabel("Geographic distance between indivudals pairs")
	plt.ylabel("Proportion of genome shared by IBD segment")
	plt.title("Relation between geographic distance and IBD shared")
	plt.savefig(args.output_folder.rstrip("/") + '/Prop_vs_distGeo_' + prefix + '.png', dpi=200, bbox_inches="tight")
