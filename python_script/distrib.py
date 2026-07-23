import argparse  
import pandas as pd
import matplotlib 
import os 
if os.environ.get("DISPLAY", "") == "":
    matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np
from pathlib import Path

### Define arguments ###
# Initialize parser
parser = argparse.ArgumentParser(description="Script to obtain the length distribution of IBD segments from a csv file")

# Add arguments with a default value or not 
parser.add_argument("--input_file", "-i", type=str,required=True, help="Path for input IBD file")
parser.add_argument("--output_folder", "-o", type=str, required=True, help="Path for output folder")
parser.add_argument("--curve", "-c", action="store_true", help="Add analytical curve described in Palamara et al (2012)")
parser.add_argument("--chr_length", "-l", type=int, default = 1e8, help="Length of simulated chromosome")
parser.add_argument("--haplotypes", "-hap", type=int, default=200, help="Number of haplotypes sampled : 2 times number of individuals sampled")
parser.add_argument("--min_length", "-l_min", type=int, default=0, help="Minimum length of represented segments")
parser.add_argument("--max_length", "-l_max", type=int, default=None, help="Maximum length of represented segments")
parser.add_argument("--min_time", "-t_min", type=int, default=None, help="Minimum age of represented segments")
parser.add_argument("--max_time", "-t_max", type=int, default=None, help="Maximum age of represented segments")
args = parser.parse_args()

# Save the prefix of the input file  
input_path = Path(args.input_file)
filename = input_path.name.split(".csv")[0]
prefix = filename.split("IBD_coord_")[1]

### Import csv ###
dfIBD = pd.read_csv(args.input_file)
print(f"Loaded a csv file containing {dfIBD.Pair.count()} IBD segments")


### Filter IBD segments by minimal and/or maximal ages, if specified ### 
if args.min_time is not None:
	dfIBD = dfIBD[dfIBD["Time"] >= args.min_time]
	
if args.max_time is not None:
	dfIBD = dfIBD[dfIBD["Time"] >= args.max_time]
	

### Return the number of IBD segments recorded for each file ###
num_segments = len(dfIBD.Pair) 
print(f"Loaded two csv files containing {num_segments} IBD segments")


### Retrieve variable values ###
L = args.chr_length
N = int(args.haplotypes / 2)

### Create histogram of IBD segments length ### 
plt.figure()
res = plt.hist(dfIBD["Length"], range = (0, L), bins=40, color = '#66CC99', edgecolor = 'white', alpha = 0.7)
plt.xlabel('Length')
plt.ylabel('Number')
plt.title(f'Length distribution of simulated IBD segments', fontsize = 16, loc = 'center')
plt.suptitle(prefix, fontsize = 12)
plt.yscale('log')
plt.xlim(args.min_length, args.max_length)
plt.savefig(args.output_folder.rstrip("/") + '/distrib_' + prefix + '.png')		
print(f"Generated IBD segments length distribution plot")
    
### Add analytical curve ###
if args.curve:																				
	# Define bounds for each column of the histogram 
	limits = list(map(int,res[1]))													  	  # convert res[1] in integer list 
	m = np.add(limits[:-1],limits[1:])/2												  # calculate mid intervals by taking succesive limits  
	M = m.tolist()																		  # convert m into a list 

	# Calculate analytical curve 
	C = []																				  # create  a list of curve values 
	for i in range(len(limits)-1):														  # loop over each limit interval (except the last, as two limits are compared at each iteration)
		u = limits[i]*10**(-6)															  # convert u limits (inferior limit) into Morgans units (1 cM = 10⁶ Morgans)
		v = (limits[i+1]-1)*10**(-6)													  # convert v limits (superior limit) into Morgans units (1 cM = 10⁶ Morgans)
		proba = 50**2*(((50+N*u)**(-2))-((50+N*v)**(-2)))							  # calculate the probability of observing a segment in this interval using a formula inspired by Palamara et al. (2012)
		C = C + [proba * num_segments]												  	  # add calculated probabilty to C list, weighted by total number of IBD segments 

	# Add analytical curve #
	plt.scatter(M, C, marker = '+', color = 'red')										  
	plt.savefig(args.output_folder.rstrip("/") + '/distrib_curve_' + prefix + '.png')										      
	print(f"Generated IBD segments length distribution plot with analytical curve")
